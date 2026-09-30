R""(
/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

void kernel jump_flooding(global const float *i_prev,
                          global const float *j_prev,
                          global float       *i_next,
                          global float       *j_next,
                          global float       *edt,
                          const int           nx,
                          const int           ny,
                          const int           step)
{
  int2 g = {get_global_id(0), get_global_id(1)};
  if (g.x >= nx || g.y >= ny) return;
  int index = linear_index(g.x, g.y, nx);

  int2 best = (int2)((int)i_prev[index], (int)j_prev[index]);

  // compute current best distance (if exists)
  float best_dist = FLT_MAX;

  if (best.x >= 0)
  {
    float dx = (float)(best.x - g.x);
    float dy = (float)(best.y - g.y);
    best_dist = dx * dx + dy * dy;
  }

  // can't do better than that...
  if (best_dist != 0.f)
  {

    // check 8 neighbors at offsets (p, q) where p,q in {-step, 0, step} except
    // (0, 0)
    for (int q = -step; q <= step; q += step)
      for (int p = -step; p <= step; p += step)
      {
        if (p == 0 && q == 0) continue;

        int ip = g.x + p;
        int jq = g.y + q;

        if (ip < 0 || ip >= nx || jq < 0 || jq >= ny) continue;

        int  nidx = linear_index(ip, jq, nx);
        int2 candidate_seed = (int2)(i_prev[nidx], j_prev[nidx]);

        if (candidate_seed.x < 0) continue; // neighbor had no seed

        float dx = (float)(candidate_seed.x - g.x);
        float dy = (float)(candidate_seed.y - g.y);
        float d2 = dx * dx + dy * dy;

        if (d2 < best_dist)
        {
          best_dist = d2;
          best = candidate_seed;
        }
      }
  }

  i_next[index] = best.x;
  j_next[index] = best.y;
  edt[index] = best_dist;
}
)""
