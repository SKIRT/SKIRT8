/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MasterSlaveCommunicator.hpp"
#include "FatalError.hpp"
#include "Parallel.hpp"
#include "ProcessManager.hpp"
#include <mutex>
#include <thread>

////////////////////////////////////////////////////////////////////

namespace
{
    // flag becomes true if the main thread has been initialized
    std::once_flag _initialized;

    // will be set to the first thread that invokes setup()
    std::thread::id _mainThread;

    // sets the main thread; will be called once (for the first thread that invokes setup()
    void setMainThread()
    {
        _mainThread = std::this_thread::get_id();
    }
}

////////////////////////////////////////////////////////////////////

MasterSlaveCommunicator::MasterSlaveCommunicator()
{
}

////////////////////////////////////////////////////////////////////

MasterSlaveCommunicator::MasterSlaveCommunicator(SimulationItem* parent)
{
    parent->addChild(this);
}

////////////////////////////////////////////////////////////////////

void MasterSlaveCommunicator::setupSelfBefore()
{
    ProcessCommunicator::setupSelfBefore();

    // initialize the main thread if needed; subsequently verify that we are running in the same thread
    std::call_once(_initialized, setMainThread);
    if (std::this_thread::get_id() != _mainThread)
        throw FATALERROR("Must be invoked from the thread that initialized the first MasterSlaveCommunicator");
}

////////////////////////////////////////////////////////////////////

MasterSlaveCommunicator::~MasterSlaveCommunicator()
{
    releaseSlaves();
}

////////////////////////////////////////////////////////////////////

void MasterSlaveCommunicator::setLocalSlaveCount(int value)
{
    if (_acquired) throw FATALERROR("Slaves are already acquired");
    _factory.setMaxThreadCount(value);
}

////////////////////////////////////////////////////////////////////

int MasterSlaveCommunicator::localSlaveCount() const
{
    return _factory.maxThreadCount();
}

////////////////////////////////////////////////////////////////////

void MasterSlaveCommunicator::setMaxMessageSize(size_t value)
{
    if (_acquired) throw FATALERROR("Slaves are already acquired");
    _bufsize = value;
}

////////////////////////////////////////////////////////////////////

size_t MasterSlaveCommunicator::maxMessageSize() const
{
    return _bufsize;
}

////////////////////////////////////////////////////////////////////

int MasterSlaveCommunicator::master() const
{
    return 0;
}

int MasterSlaveCommunicator::slave() const
{
    if (isMultiProc()) return rank();
    if (_performing) return _factory.currentThreadIndex() + 1;
    return 0;
}

////////////////////////////////////////////////////////////////////

bool MasterSlaveCommunicator::isMaster() const
{
    return !isSlave();
}

////////////////////////////////////////////////////////////////////

bool MasterSlaveCommunicator::isSlave() const
{
    return _performing || (isMultiProc() && rank());
}

////////////////////////////////////////////////////////////////////

int MasterSlaveCommunicator::registerTask(Task task)
{
    if (_acquired) throw FATALERROR("Slaves are already acquired");
    _tasks.push_back(task);
    return _tasks.size()-1;
}

////////////////////////////////////////////////////////////////////

int MasterSlaveCommunicator::taskCount() const
{
    return _tasks.size();
}

////////////////////////////////////////////////////////////////////

void MasterSlaveCommunicator::acquireSlaves()
{
    if (_acquired) throw FATALERROR("Slaves are already acquired");
    _acquired = true;
    if (isMultiProc() && isSlave())
    {
        doSlaveObeyLoop();
        _acquired = false;
    }
}

////////////////////////////////////////////////////////////////////

void MasterSlaveCommunicator::releaseSlaves()
{
    if (_performing) throw FATALERROR("Still performing tasks");
    if (isMultiProc() && _acquired && isMaster()) stopObeying();
    _acquired = false;
}

////////////////////////////////////////////////////////////////////

namespace
{
    // simple class to set a boolean flag just for the scope of the SetFlag object
    // (this is safe even when an exception is thrown)
    class SetFlag
    {
    public:
        SetFlag(bool* flag) : _flag(flag) { *_flag = true; }
        ~SetFlag() { *_flag = false; }
    private:
        bool* _flag;
    };
}

////////////////////////////////////////////////////////////////////

namespace
{
    // simple class to serve as a target for local parallel execution
    class LocalTarget : public ParallelTarget
    {
    public:
        LocalTarget(MasterSlaveCommunicator::Task& task, const vector<SerializedData>& inputVector)
            : _task(task), _inputVector(inputVector), _outputVector(inputVector.size()) { }
        void body(size_t index) { _outputVector[index] = _task(_inputVector[index]); }
        const vector<SerializedData>& outputVector() { return _outputVector; }
    private:
        MasterSlaveCommunicator::Task& _task;
        const vector<SerializedData>& _inputVector;
        vector<SerializedData> _outputVector;
    };
}

////////////////////////////////////////////////////////////////////

vector<SerializedData> MasterSlaveCommunicator::performTask(int taskIndex, const vector<SerializedData>& inputVector)
{
    if (std::this_thread::get_id() != _mainThread)
        throw FATALERROR("Must be invoked from the thread that initialized the first MasterSlaveCommunicator");
    if (_performing) throw FATALERROR("Already performing tasks");
    if (isSlave()) throw FATALERROR("Only the master can command the slaves");
    if (taskIndex < 0 || static_cast<size_t>(taskIndex) >= _tasks.size()) throw FATALERROR("Task index out of range");

    // bracket performing tasks with flag to control return value of isMaster() / isSlave()
    SetFlag flag(&_performing);

    if (isMultiProc())
    {
         return doMasterCommandLoop(taskIndex, inputVector);
    }
    else
    {
        LocalTarget target(_tasks[taskIndex], inputVector);
        _factory.parallel()->call(&target, inputVector.size());
        return target.outputVector();
    }
}

////////////////////////////////////////////////////////////////////

vector<SerializedData> MasterSlaveCommunicator::performTask(const vector<SerializedData>& inputVector)
{
    return performTask(0, inputVector);
}

////////////////////////////////////////////////////////////////////

vector<SerializedData> MasterSlaveCommunicator::doMasterCommandLoop(int taskIndex,
                                                                    const vector<SerializedData>& inputVector)
{
    // prepare an output vector of the appropriate size with data buffers of the appropriate size
    int numitems = inputVector.size();
    vector<SerializedData> outputVector(numitems);

    // prepare a vector to remember the index of most recent item handed out to each slave
    vector<int> itemForSlave(size());

    // the index of the next item to be handed out
    int numsent = 0;

    // hand out an item to each slave (unless there are less items than slaves)
    for (int slave=1; slave<size() && numsent<numitems; slave++,numsent++)
    {
        ProcessManager::sendDoubleBuffer(inputVector[numsent].data(), inputVector[numsent].used(), slave, taskIndex);
        itemForSlave[slave] = numsent;
    }

    // receive results, handing out more items until all have been handed out
    SerializedData output(_bufsize);
    for (int i=0; i<numitems; ++i)
    {
        // receive a message from any slave
        int slave;
        ProcessManager::receiveDoubleBuffer(output.data(), output.allocated(), slave);

        // put the result in the output vector
        outputVector[itemForSlave[slave]] = output;

        // if more items are available, hand one to this slave
        if (numsent<numitems)
        {
            ProcessManager::sendDoubleBuffer(inputVector[numsent].data(), inputVector[numsent].used(), slave, taskIndex);
            itemForSlave[slave] = numsent;
            numsent++;
        }
    }
    return outputVector;
}

////////////////////////////////////////////////////////////////////

void MasterSlaveCommunicator::doSlaveObeyLoop()
{
    SerializedData input(_bufsize);
    while (true)
    {
        // receive the next message from the master
        int tag;
        ProcessManager::receiveDoubleBuffer(input.data(), input.allocated(), master(), tag);

        // if the message tag specifies a non-existing task, terminate the obey loop
        if (tag < 0 || static_cast<size_t>(tag) >= _tasks.size()) break;

        // perform the requested task, deserializing and serializing QVariant from/to buffer
        SerializedData output = _tasks[tag](input);

        // send the result back to the master
        ProcessManager::sendDoubleBuffer(output.data(), output.used(), master(), tag);
    }
}

////////////////////////////////////////////////////////////////////

void MasterSlaveCommunicator::stopObeying()
{
    for (int slave=1; slave<size(); slave++)
    {
        // send an empty message with a tag that specifies a non-existing task
        SerializedData empty;
        ProcessManager::sendDoubleBuffer(empty.data(), empty.used(), slave, _tasks.size());
    }
}

////////////////////////////////////////////////////////////////////
