#include <thread>
#include <random>
#include <hrtree/isfc/key_gen.hpp>
#include <hrtree/sorting/radix_sort.hpp>
#include <hrtree/zip/zip.hpp>
#include <hrtree/zip/index_iterator.hpp>
#include <utility/rndutils.hpp>
#include "relaxation.hpp"


// TL;DR

namespace relaxation {


  // thread local, used for wiggle distribution
  static thread_local auto reng = rndutils::make_random_engine_low_entropy<>();


  int Relaxation::do_relaxation()
  {
    int k = 0;
    do {
      ++k; 
      hilbert_sort();
      build_hrtree();
    } while (query_hrtree());
    return k;
  }
  
  
  void Relaxation::build_hrtree()
  {
    rtree_.build(particles_.cbegin(), particles_.cend(), [](const auto& proxy) {
      return mbr_t(proxy.pos, proxy.radius);
    });
  }
  
  
  // sort particles by their Hilbert value
  void Relaxation::hilbert_sort()
  {
    using hrtree::zip::make_zip;
    const auto N = particles_.size();
    hrtree::key_gen<key_t, glm::vec2> kgen(area_.p0() - 2.0f, area_.p1() + 2.0f);
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
      key_[i].first = kgen(particles_[i].pos).asWord();
    }
    bool swaped = hrtree::radix_sort(make_zip(key_.begin(), particles_.begin()),
                                     make_zip(key_.end(), particles_.end()),
                                     make_zip(kbuf_.begin(), pbuf_.begin()));
    if (swaped) {
      particles_.swap(pbuf_);
      key_.swap(kbuf_);
    }
  }


  // does the heavy lifting
  // minimizing the forces between particles
  bool Relaxation::query_hrtree()
  {
    const auto N = static_cast<int>(particles_.size());
    const auto nt = (param_t::num_threads) ? param_t::num_threads : std::thread::hardware_concurrency();
    bool dirty = false;
#pragma omp parallel for schedule(static) reduction(|: dirty)
    for (int i = 0; i < N; ++i) {
      collider_t collider(static_cast<size_t>(i), particles_);
      auto mbr = mbr_t{ particles_[i].pos, particles_[i].radius - param_t::ceps };
      rtree_.query(hrtree::mbr_intersect_policy<mbr_t>(mbr), collider);
      pbuf_[i] = particles_[i];
      if (collider.dirty) {
        auto force = collider.force;
        if (glm::vec2(0) == force) {
          // we got stuck, shouldn't happen though
          auto dwiggle = std::uniform_real_distribution<float>(-param_t::wiggle, +param_t::wiggle);
          force.x = dwiggle(reng);
          force.y = dwiggle(reng);
          pbuf_[i].tex = -10.0f;
        }
        // integration step, sort of
        pbuf_[i].pos += force;
        pbuf_[i].tex += 0.002f;
        dirty = true;
      }
    }
    pbuf_.swap(particles_);
    return dirty;
  }

}
