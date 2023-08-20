#pragma once

#include <cstdint>

namespace MDPAT
{
    struct StepRange
    {
        StepRange() :
            initStep(0UL), endStep(0UL), dumpStep(0UL), nSteps(0UL)
        {}
        StepRange(uint64_t iStep, uint64_t eStep, uint64_t dStep) :
            initStep(iStep), endStep(eStep), dumpStep(dStep)
        {
            nSteps = (eStep - iStep) / dStep + 1UL;
        }
        
        uint64_t initStep = 0UL;
        uint64_t endStep = 0UL;
        uint64_t dumpStep = 0UL;
        uint64_t nSteps = 0UL;
    };

    template<typename Range>
    class StepIterator
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = uint64_t;
        using pointer = size_t;
        using reference = value_type&;
    public:
        StepIterator (pointer idx, value_type iStep, value_type eStep, value_type dStep, value_type nStep)
            : index(idx), initStep(iStep), endStep(eStep), dumpStep(dStep), nSteps(nStep), step(iStep + dStep * idx) {}

        reference operator*() const { return step; }
        // pointer operator->() { return m_ptr; }

        StepIterator& operator++()
        {
            ++index;
            step += dumpStep;
            return *this;
        }
        StepIterator operator++(int)
        {
            StepIterator tmp = *this;
            ++index;
            step += dumpStep;
            return tmp;
        }
        StepIterator& operator--()
        {
            --index;
            step -= dumpStep;
            return *this;
        }
        StepIterator operator--(int)
        {
            StepIterator tmp = *this;
            --index;
            step -= dumpStep;
            return tmp;
        }

        reference operator[](pointer idx) { return step + dumpStep * idx; }

        friend bool operator==(const Iterator& a, const Iterator& b) const { return (a.step == b.step) && (a.index == b.index); }
        friend bool operator!=(const Iterator& a, const Iterator& b) const { return (a.step != b.step) && (a.index != b.index); }
    private:
        size_t index;
        uint64_t step;
        uint64_t initStep;
        uint64_t endStep;
        uint64_t dumpStep;
        uint64_t nSteps;
    }

    class StepRangeIt
    {
    public:
        using Iterator = StepIterator<StepRangeIt>;
    public:
        StepRange() :
            initStep(0UL), endStep(0UL), dumpStep(0UL), nSteps(0UL) {}

        StepRange(uint64_t iStep, uint64_t eStep, uint64_t dStep) :
            initStep(iStep), endStep(eStep), dumpStep(dStep)
        {
            nSteps = (eStep - iStep) / dStep + 1UL;
        }

        StepIterator begin()
        {
            return Iterator(0UL, initStep, endStep, dumpStep, nSteps);
        }

        StepIterator end()
        {
            return Iterator(nSteps, initStep, endStep, dumpStep, nSteps);
        }
        
        uint64_t initStep = 0UL;
        uint64_t endStep = 0UL;
        uint64_t dumpStep = 0UL;
        uint64_t nSteps = 0UL;
    }
}
