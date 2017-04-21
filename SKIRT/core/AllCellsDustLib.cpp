/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "AllCellsDustLib.hpp"
#include "DustSystem.hpp"

////////////////////////////////////////////////////////////////////

int AllCellsDustLib::numEntries() const
{
    return find<DustSystem>()->numCells();
}

////////////////////////////////////////////////////////////////////

vector<int> AllCellsDustLib::mapping() const
{
    // create a mapping from each cell index to identical library index
    int Ncells = numEntries();
    vector<int> nv(Ncells);
    for (int m=0; m<Ncells; m++) nv[m] = m;
    return nv;
}
