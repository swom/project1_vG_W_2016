#ifndef CELLS_RELAXATION_HPP_INCLUDED
#define CELLS_RELAXATION_HPP_INCLUDED

#include <vector>
#include <hrtree/mbr_build_policy.hpp>
#include <hrtree/mbr_intersect_policy.hpp>
#include <hrtree/adapt_mbr.hpp>
#include <hrtree/isfc/hilbert.hpp>
#include <hrtree/rtree.hpp>
#include <glm/glm.hpp>
#include <glmutils/bbox.hpp>


namespace relaxation {
  using mbr_t = glmutils::bbox2;
}

// announce types we want to use to the hrtree library
HRTREE_ADAPT_POINT_FUNCTION(glm::vec2, float, 2, (float*)std::addressof);
HRTREE_ADAPT_MBR_MEMBERS(relaxation::mbr_t, glm::vec2, p0(), p1());


namespace relaxation {

  // some parameters for the relaxation step
  struct param_t {
    static constexpr float eps = 10e-6f;
    static constexpr float ceps = 0.001f;         // not required but doesn't harm
    static constexpr float spring = 1.0f;         // spring force between particles
    static constexpr float wiggle = 0.0001f;      // wiggle step if particle is trapped
    static constexpr unsigned num_threads = 0;    // 0:  std::thread::hardware_concurrency();
  };


  // proxy object 
  struct particle_t {
    glm::vec2 pos;
    float radius;
    float tex = 0.f;   // texture coordinate into color map (debugging/visualization)
  };


  // the 'QueryFun' function object used in rtree::query.
  struct collider_t
  {
    collider_t(size_t ii, const std::vector<particle_t>& particles) noexcept
      : i(ii), pp(particles.data())
    {}

    // this is called during the traversal of the Hilbert Rtree for each 
    // particle[j] whose bounding box overlaps with the bounding box of
    // particle[i]. Actually, the RTree knows about bounding boxes only,
    // thus we have to live with the index j.
    void operator()(size_t j) noexcept
    {
      if (i == j) return;
      auto& pi = pp[i];
      auto& pj = pp[j];
      const auto r_ij = pj.pos - pi.pos;
      const auto collide_dist = pi.radius + pj.radius;
      const auto dist2 = glm::dot(r_ij, r_ij);
      if (dist2 < collide_dist * collide_dist) {
        const auto rdist = 1.0f / std::sqrt(dist2 + param_t::eps * collide_dist);
        const auto n = r_ij * rdist;
        // spring force along normal ~ penetration distance
        force += (param_t::spring + param_t::ceps) * (dist2 * rdist - collide_dist) * n;
        dirty = true;
      }
    }

  public:
    const size_t i;
    const particle_t* pp;
    glm::vec2 force = glm::vec2(0);
    bool dirty = false;
  };


  class Relaxation
  {
  public:
    Relaxation(param_t param) noexcept
      : param_(param)
    {}

    // minimizes spring forces for objects of type Type in the range [first, last) 
    // 
    // ctp shall convert the input type to particle_t:
    //    particle_t c2p(const Type& a)
    // 
    // cfp is applied after relaxation:
    //    void cfp(Type& a, const particle_t& p)
    //
    // returns number of relaxation steps
    template <typename RaIt, typename CTP, typename CFP>
    int operator()(RaIt first, RaIt last, CTP&& ctp, CFP&& cfp);

  private:
    typedef hrtree::hilbert<2, 15>::type key_t;
    typedef std::pair<key_t::word_type, unsigned> keyidx_t;

    int do_relaxation();
    void build_hrtree();
    void hilbert_sort();
    bool query_hrtree();

    std::vector<particle_t> particles_, pbuf_;   // ping pong
    std::vector<keyidx_t> key_, kbuf_;           // ping pong
    hrtree::rtree<mbr_t> rtree_;
    glmutils::bbox2 area_;
    const param_t param_;
  };


  template <typename RaIt, typename CTP, typename CFP>
  int Relaxation::operator()(RaIt first, RaIt last, CTP&& ctp, CFP&& cfp)
  {
    const auto N = static_cast<uint32_t>(std::distance(first, last));
    if (N < 2) return 0;
    glmutils::bbox2 area(glm::vec2{first->get_x(),first->get_y()}, 0.f);
    particles_.resize(N); pbuf_.resize(N);
    key_.resize(N); kbuf_.resize(N);
    for (uint32_t i = 0; i < N; ++i) {
      particle_t p = ctp(*(first + i));
      key_[i].second = i;
      particles_[i] = p;
      glmutils::include(area, p.pos);
    }
    int k = do_relaxation();
    // copy back while restoring  original order
    for (uint32_t i = 0; i < N; ++i) {
      cfp(*(first + key_[i].second), particles_[i]);
    }
    return k;
  }

}

#endif

