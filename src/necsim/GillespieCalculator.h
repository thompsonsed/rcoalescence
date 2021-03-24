//
// Created by sam on 07/09/19.
//


#ifndef NECSIM_GILLESPIECALCULATOR_H
#define NECSIM_GILLESPIECALCULATOR_H

#include <map>
#include <numeric>
#include <memory>

#include "Cell.h"
#include "RNGController.h"
#include "MapLocation.h"

using std::shared_ptr;
using std::vector;

namespace necsim
{
    using namespace random_numbers;

    /**
     * @brief Container for the different event types that can occur during the Gillespie Algorithm.
     */
    enum class EventType
    {
        undefined, cell_event, map_event, sample_event
    };

    enum class CellEventType
    {
        undefined, dispersal_event, coalescence_event, speciation_event
    };

    class GillespieProbability
    {
    protected:
        double dispersal_outside_cell_probability;
        double coalescence_probability;
        double speciation_probability;
        double random_number;
        MapLocation location;

    public:

        GillespieProbability() : GillespieProbability(MapLocation())
        { }

        explicit GillespieProbability(const MapLocation &c) : dispersal_outside_cell_probability(0.0),
                                                              coalescence_probability(0.0), speciation_probability(0.0),
                                                              random_number(0.0), location(c)
        {

        }

        void setDispersalOutsideCellProbability(const double &d);

        void setCoalescenceProbability(const double &c);

        void setSpeciationProbability(const double &s);

        void setRandomNumber(const double &r);

        double getInCellProbability() const;

        CellEventType generateRandomEvent(const shared_ptr<RNGController> &rng) const;

        MapLocation &getMapLocation();

        const MapLocation &getMapLocation() const;

        /**
         * @brief Gets the parameter for the exponential function
         *
         * The rate is per birth-death event on the whole landscape.
         *
         * @param summed_death_rate
         * @param n
         * @param total_n
         * @return lambda
         */
        double getLambda(const double &local_death_rate, const double &summed_death_rate, const unsigned long &n) const;

        double calcTimeToNextEvent(const double &local_death_rate,
                                   const double &mean_death_rate,
                                   const unsigned long &n) const;

        void reset();

        friend std::ostream &operator<<(std::ostream &os, const GillespieProbability &gp);

        friend std::istream &operator>>(std::istream &is, GillespieProbability &gp);
    };

    class GillespieHeapNode
    {
    private:
        void updateHeapPosition()
        {
            // Enclosing heap vector not known or no locator to report to
            if(heap == nullptr || locator == nullptr)
            {
                return;
            }

            // Node is outside its enclosing heap, i.e. it was probably moved to a temporary place
            if(this < heap->data() || this >= (heap->data() + heap->size()))
            {
                return;
            }

            // Store the heap vector index with the locator
            *locator = this - heap->data();

        }

    public:
        Cell cell;
        double time_of_event;
        EventType event_type;

        // Pointer to index in heap
        vector<GillespieHeapNode>* heap;
        unsigned long* locator;

        GillespieHeapNode(const Cell cell,
                          const double time_of_event,
                          const EventType &e,
                          vector<GillespieHeapNode>* heap,
                          unsigned long* locator) : cell(cell), time_of_event(time_of_event), event_type(e), heap(heap),
                                                    locator(locator)
        {
            updateHeapPosition();
        }

        GillespieHeapNode(const double time_of_event, const EventType &e) : GillespieHeapNode(Cell(),
                                                                                              time_of_event,
                                                                                              e,
                                                                                              nullptr,
                                                                                              nullptr)
        { }

        GillespieHeapNode() : GillespieHeapNode(Cell(), 0.0, EventType::undefined, nullptr, nullptr)
        { }

        GillespieHeapNode(const GillespieHeapNode &other) noexcept
        {
            cell = other.cell;
            time_of_event = other.time_of_event;
            event_type = other.event_type;

            heap = other.heap;
            locator = other.locator;

            updateHeapPosition();
        }

        GillespieHeapNode &operator=(const GillespieHeapNode &other)
        {
            cell = other.cell;
            time_of_event = other.time_of_event;
            event_type = other.event_type;

            heap = other.heap;
            locator = other.locator;

            updateHeapPosition();

            return *this;
        }

        GillespieHeapNode(GillespieHeapNode &&other) noexcept
        {
            cell = other.cell;
            time_of_event = other.time_of_event;
            event_type = other.event_type;

            heap = other.heap;
            locator = other.locator;

            updateHeapPosition();
        }

        GillespieHeapNode &operator=(GillespieHeapNode &&other) noexcept
        {
            cell = other.cell;
            time_of_event = other.time_of_event;
            event_type = other.event_type;

            heap = other.heap;
            locator = other.locator;

            updateHeapPosition();

            return *this;
        }

        bool operator<(const GillespieHeapNode &n) const
        {
            return time_of_event > n.time_of_event;
        }
    };

}

#endif //NECSIM_GILLESPIECALCULATOR_H
