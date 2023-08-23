#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "stepRange.hpp"

namespace MDPAT
{

class Trajectory
{
public:
    enum class Axis {NONE = 0, FRAMES = 1, ATOMS = 2, PROPS = 3};
    typedef std::array<Axis, 3> AxisOrder;
    typedef std::array<size_t, 3> Dimensions;
public:
    Trajectory();
    ~Trajectory();

    void read(const std::filesystem::path&, const MDPAT::StepRange&);
    void read(const std::filesystem::path&);
    void read(const std::vector<std::filesystem::path>&, const MDPAT::StepRange&);
    void read(const std::vector<std::filesystem::path>&);

    const bool isLoaded() const;
    const int getColumnIndex(const char*) const;
    const bool hasColumn(const char* label) const;
    const std::vector<std::string>& getColumnLabels() const;
    const Dimensions& getAxisLengths() const;
    const Dimensions& getAxisLengthsGlobal() const;
    const AxisOrder& getAxisOrder() const;
    const std::vector<uint64_t>& getSteps const;
    const double & operator[](std::size_t idx) const;
    
    void permuteDims(const AxisOrder&);
    // void selectColumns(const std::vector<std::string> &);
    void reset();

private:
    typedef std::array<uint32_t, 3> IdxMap;
    struct TempfileHeaderResults
    {
        int startPos = 0;
        AxisOrder order;
        Dimensions dims;
    };
private:
    void initMPI();
    void checkValidAxis(const AxisOrder&) const;

    // Permuting axes
    void swap(double *, double *);
    IdxMap& getIdxMap(const AxisOrder&, const AxisOrder&) const;
    void permuteDimsLocal(
        const AxisOrder&,
        const Dimensions&,
        const IdxMap&);

    // Read text dumpfile methods
    void reserve();
    uint64_t getTimestep(std::istream&) const;
    void readDumpHeader(std::istream&);
    void skipDumpHeader(std::istream&) const;
    void readDumpBody(std::istream&, const size_t);
    void skipDumpBody(std::istream&) const;
    
    // Tempfile methods (for transposing/perumting axes)
    int writeTempfileHeader(std::ostream &outstream) const;
    void writeTempfile() const;
    TempfileHeaderResults& readTempfileHeader(std::istream &instream) const;
    void readTempfile(const AxisOrder&);
private:
    std::vector<double> m_data;  // main data

    // MPI vars
    int m_me = 0;
    int m_nprocs = 1;

    // Accessible properties (see getters above)
    bool m_loaded = false;
    bool m_initialized = false;
    AxisOrder m_axisOrder = {Axis::FRAMES, Axis::ATOMS, Axis::PROPS};
    Dimensions m_axisLengths = {0, 0, 0};
    Dimensions m_axisLengthsGlobal = {0, 0, 0};
    std::vector<std::string> m_columnLabels = {};
    // std::vector<std::string> m_originalColumnLabels;
    
    // Dumpfile vars
    uint64_t m_nframes = 0UL;
    uint64_t m_natoms = 0UL;
    uint32_t m_ncols = 0U;
    std::array<double, 6> m_box = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<uint64_t> m_stepsGlobal;
    std::vector<uint64_t> m_steps;
    std::filesystem::path m_dumpfilePath;
    std::vector<std::filesystem::path> m_dumpfilePathsVec;

    // Tempfile vars
    bool m_tempfileExists = false;
    const char* m_tempfileName = "tempfile.bin";
    std::filesystem::path m_tempfilePath;
};

}
