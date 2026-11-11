// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::ThreadedEntityIterator
 */
#ifndef EWOMS_THREADED_ENTITY_ITERATOR_HH
#define EWOMS_THREADED_ENTITY_ITERATOR_HH

#include <algorithm>
#include <atomic>
#include <mutex>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Opm {

/*!
 * \brief Provides an STL-iterator like interface to iterate over the enties of a
 *        GridView in OpenMP threaded applications
 *
 * The threads take the entities in chunks of consecutive ones, each chunk
 * under one lock, and each entity of a chunk without a lock: the work per
 * entity can vary a lot, as for a flash in two-phase cells and none in
 * single-phase ones, so the entities are distributed dynamically, but one
 * lock per entity costs more than the work on many entities.
 *
 * ATTENTION: This class must be instantiated in a sequential context!
 */
template <class GridView, int codim>
class ThreadedEntityIterator
{
    using Entity = typename GridView::template Codim<codim>::Entity;
    using EntityIterator = typename GridView::template Codim<codim>::Iterator;
public:
    //! The number of consecutive entities a thread takes at a time.
    static constexpr unsigned defaultChunkSize = 16;

    explicit ThreadedEntityIterator(const GridView& gridView,
                                    const unsigned chunkSize = defaultChunkSize)
        : sequentialIt_(gridView.template begin<codim>())
        , sequentialEnd_(gridView.template end<codim>())
        , chunkSize_(std::max(chunkSize, 1u))
        , chunks_(maxThreads_(), Chunk{sequentialEnd_, 0})
    {}

    ThreadedEntityIterator(const ThreadedEntityIterator& other) = default;

    // begin iterating over the grid in parallel
    EntityIterator beginParallel()
    {
        return takeChunk_(chunks_[threadId_()]);
    }

    // returns true if the last element was reached
    bool isFinished(const EntityIterator& it) const
    { return it == sequentialEnd_; }

    // make sure that the loop over the grid is finished
    void setFinished()
    {
        std::lock_guard<std::mutex> guard(mutex_);
        sequentialIt_ = sequentialEnd_;
        finished_.store(true, std::memory_order_relaxed);
    }

    // prefix increment: goes to the next element of the calling thread's
    // chunk, or to the first one of a chunk which is not yet worked on by any
    // thread
    EntityIterator increment()
    {
        Chunk& chunk = chunks_[threadId_()];
        if (chunk.remaining > 0 && !finished_.load(std::memory_order_relaxed)) {
            --chunk.remaining;
            return ++chunk.it;
        }
        return takeChunk_(chunk);
    }

private:
    // The entity a thread works on and the number of entities after it in its chunk.
    struct alignas(64) Chunk
    {
        EntityIterator it;
        unsigned remaining;
    };

    EntityIterator takeChunk_(Chunk& chunk)
    {
        std::lock_guard<std::mutex> guard(mutex_);
        chunk.it = sequentialIt_;
        unsigned taken = 0;
        while (taken < chunkSize_ && sequentialIt_ != sequentialEnd_) {
            ++sequentialIt_;
            ++taken;
        }
        chunk.remaining = taken > 0 ? taken - 1 : 0;
        return chunk.it;
    }

    static unsigned maxThreads_()
    {
#ifdef _OPENMP
        return static_cast<unsigned>(std::max(omp_get_max_threads(), 1));
#else
        return 1;
#endif
    }

    static unsigned threadId_()
    {
#ifdef _OPENMP
        return static_cast<unsigned>(omp_get_thread_num());
#else
        return 0;
#endif
    }

    EntityIterator sequentialIt_;
    EntityIterator sequentialEnd_;
    unsigned chunkSize_;
    std::vector<Chunk> chunks_;
    std::atomic<bool> finished_{false};

    std::mutex mutex_;
};

} // namespace Opm

#endif
