#include "trajectory.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <mpi.h>

#ifndef OMPI_MPI_H
#include "mpi_stub.h"
#endif

#include "error.hpp"
#include "splitValues.hpp"
#include "mdbin.h"

using std::string;
using std::vector;

namespace MDPAT
{

void Trajectory::initMPI()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &m_me);
    MPI_Comm_size(MPI_COMM_WORLD, &m_nprocs);
}

const bool Trajectory::isLoaded() const
{
    return m_loaded;
}

const int Trajectory::getColumnIndex(const char* label) const
{
    auto it = std::find(m_columnLabels.begin(), m_columnLabels.end(), label);
    return 0;
}

const bool Trajectory::hasColumn(const char* label) const
{
    return m_columnLabels.end() != std::find(m_columnLabels.begin(), m_columnLabels.end(), label);
}

const std::vector<std::string>& Trajectory::getColumnLabels() const
{
    return m_columnLabels;
}

const Trajectory::Dimensions& Trajectory::getAxisLengths() const
{
    return m_axisLengths;
}

const Trajectory::Dimensions& Trajectory::getAxisLengthsGlobal() const
{
    return m_axisLengthsGlobal;
}

const Trajectory::AxisOrder& Trajectory::getAxisOrder() const
{
    return m_axisOrder;
}

const double & Trajectory::operator[](std::size_t idx) const
{
    return m_data[idx];
}

void Trajectory::read(
    const std::filesystem::path& dumpfile,
    const StepRange& stepRange)
{
    // Remove tempfile, if it exists
    if (m_me == 0 && std::filesystem::is_regular_file(m_tempfilePath))
        std::filesystem::remove(m_tempfilePath);

    // Calculate first frame, num frames from step range
    m_stepRange = stepRange;
    const auto [firstFrame, numFrames] = splitValues(m_stepRange.nSteps, m_me, m_nprocs);

    // Assume dumpfile exists and open, start reading, throw if it doesn't
    std::ifstream instream(dumpfile);
    if (!instream.good())
        errorAll(Error::IOERROR, "Could not open file %s", m_dumpfilePath.c_str());
    
    // Loop through my expected timesteps
    bool headerDone = false;
    uint64_t firstStep = firstFrame * m_stepRange.dumpStep + m_stepRange.initStep;
    uint64_t lastStep = firstStep + (numFrames - 1) * m_stepRange.dumpStep;
    uint64_t step = firstStep;
    while (step <= lastStep && instream.good())
    {
        uint64_t dumpfileTimestep = getTimestep(instream);
        if (dumpfileTimestep == step)
        {
            if (headerDone)
            {
                skipDumpHeader(instream);
            }
            else
            {
                readDumpHeader(instream);
                headerDone = true;
                reserve();
                m_data.resize(numFrames * m_natoms * m_ncols);
            }
            readDumpBody(instream);
            step += m_stepRange.dumpStep;
        }
        else if (dumpfileTimestep < step)
        {
            skipDumpHeader(instream);
            skipDumpBody(instream);
        }
        else
        {
            errorAll(Error::IOERROR, "Specified timestep %ul not found in dump file", step);
        }
    }

    // If all successful, set member vars
    MPI_Barrier(MPI_COMM_WORLD);
    m_loaded = true;

    m_axisOrder = {Trajectory::Axis::FRAMES, Trajectory::Axis::ATOMS, Trajectory::Axis::PROPS};
    m_axisLengths = {numFrames, m_natoms, m_columnLabels.size()};
    m_axisLengthsGlobal = {m_stepRange.nSteps, m_natoms, m_columnLabels.size()};

    instream.close();
}

void Trajectory::swap(double* x, double* y)
{
    double tmp = *x;
    *x = *y;
    *y = tmp;
}

/*
 Returns idxMap, a permutation of {0, 1, 2}, such that order1[i] == order2[idxMap[i]]
*/
Trajectory::IdxMap& Trajectory::getIdxMap(const Trajectory::AxisOrder& order1,
                                          const Trajectory::AxisOrder& order2) const
{
    Trajectory::IdxMap idxMap = {0U, 0U, 0U};
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            if (order1[i] == order2[j])
            {
                idxMap[i] = j;
                break;
            }
        }
    }
    return idxMap;
}

void Trajectory::checkValidAxis(const AxisOrder& axisOrder) const
{
    if (axisOrder.end() == std::find(axisOrder.begin(), axisOrder.end(), Axis::FRAMES))
        errorAll(Error::ARGUMENTERROR, "Argument axisOrder doesn't specify axis FRAMES.");
    if (axisOrder.end() == std::find(axisOrder.begin(), axisOrder.end(), Axis::ATOMS))
        errorAll(Error::ARGUMENTERROR, "Argument axisOrder doesn't specify axis ATOMS.");
    if (axisOrder.end() == std::find(axisOrder.begin(), axisOrder.end(), Axis::PROPS))
        errorAll(Error::ARGUMENTERROR, "Argument axisOrder doesn't specify axis PROPS.");

    return;
}

void Trajectory::reserve()
{
    uint64_t max_nSteps = (uint64_t)ceil((double)m_stepRange.nSteps / (double)m_nprocs);
    uint64_t max_nAtoms = m_nprocs * (uint64_t)ceil((double)m_natoms / (double)m_nprocs);
    uint64_t max_nCols = m_nprocs * (uint64_t)ceil((double)m_ncols / (double)m_nprocs);
    m_data.reserve(max_nSteps * max_nAtoms * max_nCols);
}

uint64_t Trajectory::getTimestep(std::istream &is) const
{
    string word;
    is >> word;
    if (word != "ITEM:")
        errorAll(Error::SYNTAXERROR, "Syntax error while reading dump file");
    is >> word;
    if (word != "TIMESTEP")
        errorAll(Error::SYNTAXERROR, "Syntax error while reading dump file");
    
    uint64_t timestep = 0UL;
    is >> timestep;
    return timestep;
}

void Trajectory::readDumpHeader(std::istream &is)
{
    string word;
    m_box.fill(0.0);
    m_natoms = 0UL;
    is >> word;

    if (word != "ITEM:")
        errorAll(Error::SYNTAXERROR, "Syntax error while reading dump file");

    while (is.good())
    {
        is >> word;
        if (word == "TIMESTEP")
        {
            is >> word;
        }
        else if (word == "BOX")
        {
            is >> word >> word >> word >> word; // ITEM: BOX BOUNDS ab ab ab
            for (size_t i = 0; i < m_box.size(); ++i)
                is >> m_box[i];
        }
        else if (word == "NUMBER")
        {
            is >> word >> word; // ITEM: NUMBER OF ATOMS
            is >> m_natoms;
        }
        else if (word == "ATOMS")
        {
            std::getline(is, word); // ITEM: ATOMS id type x y z ... => word == {"id", "type", ...}
            std::istringstream labelstream(word);
            m_columnLabels.clear();
            labelstream >> word;
            if (word != "id")
                errorAll(Error::SYNTAXERROR, "First column of dump file must be 'id'");
            while (labelstream.good() && word != "")
            {
                m_columnLabels.push_back(word);
                labelstream >> word;
            }
            m_ncols = m_columnLabels.size();
            return;
        }
        else
        {
            while (word != "ITEM:" && is.good())
                is >> word;
            if (!is.good())
                errorAll(Error::IOERROR, "File ended before the header finished");
        }
        if (word != "ITEM:")
            is >> word;
        if (word != "ITEM:")
            errorAll(Error::SYNTAXERROR, "Syntax error while reading dump file");
    }
    errorAll(Error::IOERROR, "File ended before the header finished");
}

void Trajectory::skipDumpHeader(std::istream& is) const
{
    string word;
    is >> word;
    while (is.good())
    {
        if (word != "ITEM:")
            continue;
        is >> word;
        if (word == "ATOMS")
        {
            is.ignore(512, '\n');
            return;
        }
    }
    errorAll(Error::IOERROR, "File ended before the header finished");
}

void Trajectory::readDumpBody(std::istream& is, const size_t offset)
{
    size_t id;
    size_t num_cols = m_columnLabels.size();

    for (size_t i = 0; i < m_natoms; ++i)
    {
        is >> id;
        for (size_t j = 0; j < num_cols; ++j)
            is >> m_data[offset + (id-1) * num_cols + j];
    }
}

void Trajectory::skipDumpBody(std::istream& is) const
{
    for (uint64_t i = 0; i < m_natoms; ++i)
        is.ignore(512, '\n');
}

int Trajectory::writeTempfileHeader(std::ostream& outstream) const
{
    outstream.write("P", 1);
    outstream.write(reinterpret_cast<const char*>(&m_nprocs), sizeof(int));
    outstream.write("R3A", 3);
    for (size_t i = 0; i < 3; ++i)
    {
        switch (m_axisOrder[i])
        {
        case Trajectory::Axis::FRAMES:
            outstream.write("F", 1);
            break;

        case Trajectory::Axis::ATOMS:
            outstream.write("A", 1);
            break;

        case Trajectory::Axis::PROPS:
            outstream.write("P", 1);
            break;
        
        default:
            break;
        }
    }
    outstream.write("D", 1);
    outstream.write(reinterpret_cast<const char*>(&m_axisLengthsGlobal), 3*sizeof(size_t));
    int curPos = outstream.tellp();
    return curPos;
}

void Trajectory::writeTempfile() const
{
    std::ofstream outstream;
    outstream.open(m_tempfileName, std::ios::binary);

    int curPos = 0;
    if (m_me == 0)
        curPos = writeTempfileHeader(outstream);

    for (size_t i = 0; i < m_nprocs; ++i)
    {
        if (m_me == i)
        {
            outstream.seekp(curPos);
            outstream.write(reinterpret_cast<const char*>(&m_data), m_data.size()*sizeof(double));
            curPos = outstream.tellp();
            if (i + 1 < m_nprocs)
                MPI_Send(&curPos, 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD);
        }
        else if (m_me == i+1)
        {
            MPI_Recv(&curPos, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // MPI_Bcast(&curPos, 1, MPI_INT, me, MPI_COMM_WORLD);
    }
    outstream.close();
}

/*
 Returns a struct with the following members:
 int startPos is the starting position for the data
 AxisOrder order is the order of the axes as written
 Dimensions dims is the length of each dimension
 Returns with startPos = -1 on error.
*/
Trajectory::TempfileHeaderResults& Trajectory::readTempfileHeader(std::istream& instream) const
{
    Trajectory::TempfileHeaderResults results;
    results.startPos = -1;
    char buf[16];
    instream.read(buf, 1);
    if (buf[0] != 'P')
        return results;

    instream.ignore(sizeof(int) + 3);
    for (size_t i = 0; i < 3; ++i)
    {
        char ax;
        instream.read(&ax, 1);
        switch (ax)
        {
        case 'F':
            results.order[i] = Trajectory::Axis::FRAMES;
            break;
        
        case 'A':
            results.order[i] = Trajectory::Axis::ATOMS;
            break;
        
        case 'P':
            results.order[i] = Trajectory::Axis::PROPS;
            break;
        
        default:
            break;
        }
    }

    uint64_t dims[] = {0UL, 0UL, 0UL};
    instream.read(reinterpret_cast<char*>(&dims), 3*sizeof(uint64_t));
    for (size_t i = 0; i < 3; ++i)
        results.dims[i] = dims[i];
    
    results.startPos = instream.tellg();
    return results;
}

void Trajectory::readTempfile(const Trajectory::AxisOrder& newAxisOrder)
{
    std::ifstream instream;
    instream.open(m_tempfileName, std::ios::binary);

    auto results = readTempfileHeader(instream);
    if (results.startPos == -1)
        errorAll(Error::IOERROR, "ERROR: Error reading file %s", m_tempfileName);

    // Now, we can say newAxisOrder[i] == results.order[new2oldIdx[i]]
    // and newGlobalLengths[i] == m_axisLengthsGlobal[new2oldIdx[i]]
    auto new2oldIdx = getIdxMap(newAxisOrder, results.order);
    const auto [myFirstIdx, nValues] = splitValues(results.dims[new2oldIdx[0]], m_me, m_nprocs);
    const uint64_t myLastIdx = myFirstIdx + nValues - 1;
    MPI_Barrier(MPI_COMM_WORLD);

    uint64_t linearReadIdx = 0UL;
    Trajectory::Dimensions arrayReadIdx = {0UL, 0UL, 0UL};
    for (arrayReadIdx[0] = 0; arrayReadIdx[0] < results.dims[0]; ++arrayReadIdx[0])
    {
        for (arrayReadIdx[1] = 0; arrayReadIdx[1] < results.dims[1]; ++arrayReadIdx[1])
        {
            for (arrayReadIdx[2] = 0; arrayReadIdx[2] < results.dims[2]; ++arrayReadIdx[2])
            {
                if (myFirstIdx <= arrayReadIdx[new2oldIdx[0]] && arrayReadIdx[new2oldIdx[0]] <= myLastIdx)
                {
                    Trajectory::Dimensions arrayWriteIdx = {
                        arrayReadIdx[new2oldIdx[0]],
                        arrayReadIdx[new2oldIdx[1]],
                        arrayReadIdx[new2oldIdx[2]]
                    };

                    uint64_t linearWriteIdx = arrayWriteIdx[0] * results.dims[new2oldIdx[1]] * results.dims[new2oldIdx[2]]
                                            + arrayWriteIdx[1] * results.dims[new2oldIdx[2]]
                                            + arrayWriteIdx[2];

                    instream.read(reinterpret_cast<char*>(&m_data[linearWriteIdx]), sizeof(double));
                }
                else
                {
                    instream.ignore(sizeof(double));
                }
            }
        }
    }
}

void Trajectory::permuteDimsLocal(
    const Trajectory::AxisOrder& newAxisOrder,
    const Trajectory::Dimensions& newLengths,
    const Trajectory::IdxMap& old2newIdx)
{
    int oldArrayIdx[3] = {0, 0, 0};
    int newArrayIdx[3] = {0, 0, 0};

    vector<bool> completed(m_data.size() - 1);
    completed[0] = true;

    // The first and last elements won't move during transposition
    for (int start = 1; start < m_data.size() - 1; ++start)
    {
        if (completed[start])
            continue;

        int newLinearIdx = start;
        // Can likely improve this loop. std::div twice on every element is bad
        while (true)
        {
            // int q = 0;
            // [q, oldArrayIdx[2]] = std::div(newLinearIdx, m_axisLengths[2]);
            const auto& div1 = std::div(newLinearIdx, m_axisLengths[2]);
            // [oldArrayIdx[0], oldArrayIdx[1]] = std::div(q, m_axisLengths[1]);
            const auto& div2 = std::div(div1.quot, m_axisLengths[1]);

            oldArrayIdx[0] = div2.quot;
            oldArrayIdx[1] = div2.rem;
            oldArrayIdx[2] = div1.rem;

            for(int i = 0; i < 3; ++i)
                newArrayIdx[old2newIdx[i]] = oldArrayIdx[i];

            newLinearIdx = newArrayIdx[0] * newLengths[1] * newLengths[2]
                         + newArrayIdx[1] * newLengths[2]
                         + newArrayIdx[2];

            completed[newLinearIdx] = true;

            if (newLinearIdx == start)
                break;
            
            swap(&m_data[start], &m_data[newLinearIdx]);
        }        
    }
}

void Trajectory::permuteDims(const AxisOrder& newAxisOrder)
{
    checkValidAxis(newAxisOrder);
    
    if (newAxisOrder == m_axisOrder)
        return;

    auto old2newIdx = getIdxMap(m_axisOrder, newAxisOrder);
    /*
     m_axisOrder[i] == newAxisOrder[idxMap[i]]
    */

    Trajectory::Dimensions newLengths = {0, 0, 0};
    for (size_t i = 0; i < 3; ++i)
        newLengths[old2newIdx[i]] = m_axisLengths[i];
    
    if (old2newIdx[0] == 0)
    {
        permuteDimsLocal(newAxisOrder, newLengths, old2newIdx);
        for (size_t i = 0; i < 3; ++i)
        {
            m_axisLengths[i] = newLengths[i];
            m_axisOrder[i] = newAxisOrder[i];
        }
        return;
    }
    else
    {
        if (!m_tempfileExists)
        {
            MDBIN::HeaderInfo info;
            info.initStep = m_stepRange.initStep;
            info.endStep = m_stepRange.endStep;
            info.deltaStep = m_stepRange.dumpStep;
            info.columnLabels = m_columnLabels;
            info.numFrames = m_axisLengthsGlobal[0];
            info.numAtoms = m_axisLengthsGlobal[1];
            info.numCols = m_axisLengthsGlobal[2];
            if (m_me == 0 && MDBIN::write(m_tempfilePath, info, m_data))
                errorAll(Error::IOERROR, "Couldn't write to file %s", m_tempfileName);
            
            MPI_Barrier(MPI_COMM_WORLD);
            for (int i = 0; i < m_nprocs; ++i)
            {
                if (i != 0 && m_me == i)
                    MDBIN::append(m_tempfilePath, m_data);

                MPI_Barrier(MPI_COMM_WORLD);
            }
            // writeTempfile();
            m_tempfileExists = true;
        }
        readTempfile(newAxisOrder);
    }

    // Below is the original implementation plan with minimal wasted memory and no writes to disk

    // Determine initial, intermediate (maximum), and final dimensions of vector
    // Resize vector to maximum necessary size
    // Copy data to appropriate buffer locations
    // MPI All-to-All data communication
    // Copy data back to appropriate locations
    // Resize vector to final size
    // Update axisOder, axisLengthsLocal, m_axisLengthsGlobal, etc.
}

void Trajectory::reset()
{

}

}