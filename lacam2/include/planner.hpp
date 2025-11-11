/*
 * lacam-star
 */

#pragma once

#include "dist_table.hpp"
#include "graph.hpp"
#include "instance.hpp"
#include "utils.hpp"

// objective function
enum Objective { OBJ_NONE, OBJ_MAKESPAN, OBJ_SUM_OF_LOSS };
std::ostream& operator<<(std::ostream& os, const Objective objective);

// PIBT agent
struct Agent {
  const uint id;
  Vertex* v_now;   // current location
  Vertex* v_next;  // next location
  Agent(uint _id) : id(_id), v_now(nullptr), v_next(nullptr) {}
};
using Agents = std::vector<Agent*>;

// low-level node
struct LNode {
  std::vector<uint> who;
  Vertices where;
  const uint depth;
  LNode(LNode* parent = nullptr, uint i = 0,
        Vertex* v = nullptr);  // who and where
};

// high-level node
struct HNode {
  static uint HNODE_CNT;  // count #(high-level node)
  const Config C;

  // tree
  HNode* parent;
  std::set<HNode*> neighbor;

  // costs
  uint g;        // g-value (might be updated)
  const uint h;  // h-value
  uint f;        // g + h (might be updated)

  // for low-level search
  std::vector<float> priorities;
  std::vector<uint> order;
  std::queue<LNode*> search_tree;

  HNode(const Config& _C, DistTable& D, HNode* _parent, const uint _g,
        const uint _h);
  ~HNode();
};
using HNodes = std::vector<HNode*>;


struct Group {
  // Parallel arrays: ids[k] has location locs[k].
  std::vector<uint> ids;     // sorted, unique
  std::vector<uint> locs;    // aligned with ids

  // Timesteps: caller appends strictly increasing.
  std::vector<uint> timesteps;

  // Optional: handle to your low-level tree if you keep it here.
  std::queue<LNode*> search_tree;

  Group() = default;

  // Construct from (id,loc) pairs in ANY order.
  Group(std::vector<std::pair<uint,uint>> pairs):
  search_tree(std::queue<LNode*>())
  {
    if (!pairs.empty()) 
    {
      // sort by id
      std::sort(pairs.begin(), pairs.end(),
                [](auto& a, auto& b){ return a.first < b.first; });
      // split
      ids.resize(pairs.size());
      locs.resize(pairs.size());
      for (size_t i = 0; i < pairs.size(); ++i) {
        ids[i]  = pairs[i].first;
        locs[i] = pairs[i].second;
      }
    }
    search_tree.push(new LNode());
  }

  // Append strictly increasing timestep (asserts in debug).
  void add_timestep(uint t) 
  {
    timesteps.push_back(t);
  }

  // Equality: same ids and aligned locations.
  bool equals(const Group& other) const 
  {
    if (ids.size() != other.ids.size()) 
      return false;
    if (locs.size() != other.locs.size()) 
      return false;
    for (size_t i = 0; i < ids.size(); i++) 
    {
      if (ids[i] != other.ids[i]) 
        return false;
      if (locs[i] != other.locs[i]) 
        return false;
    }
    return true;
  }
  bool operator==(const Group& rhs) const { return equals(rhs); }
  bool operator!=(const Group& rhs) const { return !equals(rhs); }

  // Utilities
  size_t size() const { return ids.size(); }
  bool empty()  const { return ids.empty(); }
  int   min_id() const { return empty() ? -1 : static_cast<int>(ids.front()); }

  // Reserve if you know approximate size (saves reallocations).
  void reserve(size_t n) { ids.reserve(n); locs.reserve(n); }
};

// robust 64-bit mixer
static inline uint64_t mix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

struct GroupHash {
  size_t operator()(const Group& g) const noexcept {
    uint64_t h = 0x1234567890abcdefULL;
    for (uint v : g.ids)  h ^= mix64(v) + (h<<6) + (h>>2);
    h ^= 0xf00df00df00df00dULL; // delimiter between ids and locs
    for (uint v : g.locs) h ^= mix64(v) + (h<<6) + (h>>2);
    return static_cast<size_t>(h);
  }
};

struct Planner {
  const Instance* ins;
  const Deadline* deadline;
  std::mt19937* MT;
  const int verbose;

  // hyper parameters
  const Objective objective;
  const float RESTART_RATE;  // random restart

  // solver utils
  const uint N;       // number of agents
  const uint V_size;  // number o vertices
  DistTable D;
  uint loop_cnt;      // auxiliary

  // used in PIBT
  std::vector<std::array<Vertex*, 5> > C_next;  // next locations, used in PIBT
  std::vector<float> tie_breakers;              // random values, used in PIBT
  Agents A;
  Agents occupied_now;                          // for quick collision checking
  Agents occupied_next;                         // for quick collision checking

  Planner(const Instance* _ins, const Deadline* _deadline, std::mt19937* _MT,
          const int _verbose = 0,
          // other parameters
          const Objective _objective = OBJ_NONE,
          const float _restart_rate = 0.001);
  ~Planner();
  Solution solve(std::string& additional_info);
  Solution solve_group_pibt(std::string& additional_info);
  void expand_lowlevel_tree(HNode* H, LNode* L);
  void rewrite(HNode* H_from, HNode* T, HNode* H_goal,
               std::stack<HNode*>& OPEN);
  uint get_edge_cost(const Config& C1, const Config& C2);
  uint get_edge_cost(HNode* H_from, HNode* H_to);
  uint get_h_value(const Config& C);
  bool get_new_config(HNode* H, LNode* L);
  bool funcPIBT(Agent* ai);

  // swap operation
  Agent* swap_possible_and_required(Agent* ai);
  bool is_swap_required(const uint pusher, const uint puller,
                        Vertex* v_pusher_origin, Vertex* v_puller_origin);
  bool is_swap_possible(Vertex* v_pusher_origin, Vertex* v_puller_origin);

  // utilities
  template <typename... Body>
  void solver_info(const int level, Body&&... body)
  {
    if (verbose < level) return;
    std::cout << "elapsed:" << std::setw(6) << elapsed_ms(deadline) << "ms"
              << "  loop_cnt:" << std::setw(8) << loop_cnt
              << "  node_cnt:" << std::setw(8) << HNode::HNODE_CNT << "\t";
    info(level, verbose, (body)...);
  }
};
