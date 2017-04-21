/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef MASTERSLAVECOMMUNICATOR_HPP
#define MASTERSLAVECOMMUNICATOR_HPP

#include "ParallelFactory.hpp"
#include "ProcessCommunicator.hpp"
#include "SerializedData.hpp"
#include <functional>

////////////////////////////////////////////////////////////////////

/**
An object of the MasterSlaveCommunicator class, inheriting from ProcessCommunicator, represents
an environment or ensemble of processes, which are able to communicate according to the
master-slave model. This means that all communications go through a single process, called the
master. This process sends messages to and receives messages from the other processes - the
slaves. The MasterSlaveCommunicator obtains its process environment with the setup of its base
class. Either the number of processes is one (singleprocessing mode) or is greater than one
(multiprocessing mode). In multiprocessing mode, the parallel tasks are handed out by the
master between its slave processes. In singleprocessing mode, the MasterSlaveCommunicator uses
multithreading for its parallel tasks.

<B>Usage example</B>

Here is a complete example of using this class:

\code
class Compute
{
public:
    Compute(int size, double factor)
        : _size(size), _factor(factor)
    {
        _comm.setLocalSlaveCount(4);
        _comm.registerTask([this] (const SerializedData& input) { return doSingle(input); });
    }
    void setup()
    {
        _comm.setup();
        _comm.acquireSlaves();
    }
    ~Compute()
    {
        _comm.releaseSlaves();
    }
    void doIt()
    {
        if (_comm.isMaster())
        {
            vector<SerializedData> datav(_size);
            for (int i=0; i<_size; i++) datav[i].push(static_cast<double>(i));
            datav = _comm.performTask(datav);
            for (int i=0; i<_size; i++) std::cout << datav[i].pop() << " ";
            std::cout << std::endl;
        }
    }
private:
    SerializedData doSingle(const SerializedData& input)
    {
        SerializedData data = input;
        data.push(data.pop()*_factor);
        return data;
    }
    int _size;
    double _factor;
    MasterSlaveCommunicator _comm;
};

int main(int argc, char** argv)
{
    ProcessManager pm(&argc, &argv);
    System system(argc, argv);

    {
        Compute c(7, 2.);
        c.setup();
        c.doIt();
    }
    {
        Compute c(17, 0.5);
        c.setup();
        c.doIt();
    }
}
\endcode

When executed this example will print:

\verbatim
0 2 4 6 8 10 12
0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8
\endverbatim

The constructor of the Compute class initializes a MasterSlaveCommunicator instance and
registers one of its own methods as a task to be performed by one or more slaves in parallel.
The setup() function prepares the slaves to receive work. This is split off from the
constructor to ensure that the Compute object has been fully constructed before calling the
acquireSlaves() function, which doesn't return immediately on slave processes. The doIt()
function does work only if it is running from the master; in that case it prepares some data,
commands the slaves to perform parallel work on each of the data elements, and finally prints
the results.

The main() function constructs two different Compute objects and performs the corresponding
calculations. Any context created before calling the acquireSlaves() function (such as the
multiplication factor in the example) can be used in the slave computations, even if they are
performed by multiple processes. This is possible because each slave process also executes all
code until the acquireSlaves() function is called. For a slave process, this function starts a
communication loop that exits when the master calls the releaseSlaves() function. In fact, the
two code segments are each surrounded by curly braces so that the Compute object would go out
of scope and execute its destructor, releasing the slaves.

<B>Parallel modes</B>

In singleprocessing mode, the class uses a Parallel object spawned from a privately owned
ParallelFactory instance. In multiprocessing mode, the implementation relies on MPI (Message
Passing Interface) for communication between master and slaves. Whether singleprocessing or
multiprocessing mode is used, is determined during the setup of the ProcessManager base class.
Only if MPI is available and the number of MPI processes is greater than one, multiprocessing
mode is used. In all other cases local mode is used instead.

<B>Passing data to and from the slaves</B>

Except for any context established before calling the acquireSlaves() function, all information
passed from the master to a slave and back must be serialized to a flat byte stream (at least
in multiprocessing mode). This is accomplished through the use of SerializedData objects, which
can hold an arbitrary sequence of data items of type double or vector<double>, and which can be
easily transmitted over MPI.

<B>Thread safety (or lack thereof)</B>

With the exception of isMaster() and isSlave(), all MasterSlaveCommunicator functions
(including instance construction and setup) must be invoked from the very same thread. Usually
this will be the main program thread. In some key places, a fatal error is thrown if this
restriction is violated.
*/
class MasterSlaveCommunicator : public ProcessCommunicator
{
    //============= Construction - Setup - Destruction =============

public:
    /** The default constructor creates a master-slave communicator that is \em not hooked up in a
        simulation hierarchy. The setup() function is \em not called by this constructor. */
    explicit MasterSlaveCommunicator();

    /** This constructor creates a master-slave communicator that is hooked up as a child to the
        specified parent in the simulation hierarchy, so that it will automatically be deleted. The
        setup() function is \em not called by this constructor. */
    explicit MasterSlaveCommunicator(SimulationItem* parent);

protected:
    /** This function checks whether it is invoked from the main thread. */
    void setupSelfBefore() override;

public:
    /** Releases the slaves, if applicable, and destructs the object. */
    ~MasterSlaveCommunicator();

    //====================== Other Functions =======================

public:
    /** Sets the number of slaves to be used when operating in local mode; this number is ignored
        when operating in multiprocessing mode. Throws a fatal error if called while slaves are
        acquired. */
    void setLocalSlaveCount(int value);

    /** Returns the number of slaves to be used when operating in local mode. */
    int localSlaveCount() const;

    /** Sets the maximum size of a message (as a number of double values) exchanged between master
        and slave when operating in multiprocessing mode. This number is ignored when operating in
        singleprocessing mode. The number must be large enough to accomodate any of the
        SerializedData objects passed to or returned from the registerTask() function. The default
        value is 512 doubles, which is sufficient in most cases. */
    void setMaxMessageSize(size_t value);

    /** Returns the maximum size of a message (as a number of double values) exchanged between
        master and slave when operating in multiprocessing mode. */
    size_t maxMessageSize() const;

    /** Returns the rank of the master process. */
    int master() const;

    /** Returns the rank of the calling slave thread or process (between 1 and the number of slaves),
        or zero if the caller is the master. */
    int slave() const;

    /** Returns true if the caller is the master. */
    bool isMaster() const;

    /** Returns true if the caller is a slave. */
    bool isSlave() const;

    /** Definition of the type of a function representing a task, taking serialized data as input
        and generating serialized data as output. */
    using Task = std::function<SerializedData (const SerializedData& input)>;

    /** Registers the specified function as a task. Task indices are assigned in increasing order
        starting from zero. Throws a fatal error if called while slaves are acquired. */
    int registerTask(Task task);

    /** Returns number of tasks. */
    int taskCount() const;

    /** Ensures that master and slaves are ready to command and perform tasks. In multiprocessing
        mode, the slaves block until releaseSlaves() is called. Throws a fatal error if called
        while slaves are already acquired. */
    void acquireSlaves();

    /** Release the slaves, if applicable. Does nothing if the slaves are not acquired, or if
        called from a slave. Throws a fatal error if tasks are still being performed. */
    void releaseSlaves();

    /** Make the slaves perform the task with specified index on each of the data items in the
        specified vector (in parallel). The results are returned in a vector with the same size as
        the input vector. Throws a fatal error if called while slaves are not acquired, if called
        from a slave, or if the task index is out of range. */
    vector<SerializedData> performTask(int taskIndex, const vector<SerializedData>& inputVector);

    /** Make the slaves perform the task with index zero on each of the data items in the specified
        vector. Invokes the general performTask() function with a task index of zero. */
    vector<SerializedData> performTask(const vector<SerializedData>& data);

    //====== Private Functions for multiprocessing Operation =======

private:
    /** Implements the command loop for the master process. */
    vector<SerializedData> doMasterCommandLoop(int taskIndex, const vector<SerializedData>& inputVector);

    /** Implements the obey loop for a slave process. */
    void doSlaveObeyLoop();

    /** Makes the slave processes exit their obey loop. */
    void stopObeying();

    //======================== Data Members ========================

private:
    bool _acquired{false};      // true if slaves are acquired
    bool _performing{false};    // true if we're performing a set of tasks
    ParallelFactory _factory;   // the factory used to spawn objects for local parallellization
    vector<Task> _tasks;        // registered tasks, in index order
    size_t _bufsize{512};       // the maximum message size, in doubles (for multiprocessing mode)
};

////////////////////////////////////////////////////////////////////

#endif
