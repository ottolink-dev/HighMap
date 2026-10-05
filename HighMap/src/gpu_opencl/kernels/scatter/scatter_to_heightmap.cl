R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

void kernel scatter_to_heightmap(global float       *heightmap,
                                 global float       *rock_map,
                                 global const float *items_pos_rad,
                                 global const uint  *items_class,
                                 const int           num_items,
                                 const int           nx,
                                 const int           ny,
                                 const float4        bbox,
                                 const int           shape_type,
                                 const float         height_radius_ratio,
                                 const uint          seed,
                                 const int           filter_class,
                                 const int           has_rock_map)
{
  int item_id = get_global_id(0);
  if (item_id >= num_items) return;

  uint item_class = items_class[item_id];
  if (filter_class >= 0 && item_class != (uint)filter_class) return;

  float cx = items_pos_rad[item_id * 4 + 0];
  float cy = items_pos_rad[item_id * 4 + 1];
  float cz = items_pos_rad[item_id * 4 + 2];
  float radius = items_pos_rad[item_id * 4 + 3];

  if (radius <= 0.0f) return;

  float xmin = bbox.x;
  float xmax = bbox.y;
  float ymin = bbox.z;
  float ymax = bbox.w;

  float dx_dom = xmax - xmin;
  float dy_dom = ymax - ymin;
  if (dx_dom <= 0.0f || dy_dom <= 0.0f) return;

  float px = dx_dom / (float)nx;
  float py = dy_dom / (float)ny;

  // determine bounding box in pixel grid
  int ix_min = clamp((int)floor((cx - radius - xmin) / px), 0, nx - 1);
  int ix_max = clamp((int)ceil((cx + radius - xmin) / px), 0, nx - 1);
  int iy_min = clamp((int)floor((cy - radius - ymin) / py), 0, ny - 1);
  int iy_max = clamp((int)ceil((cy + radius - ymin) / py), 0, ny - 1);

  float h_max = radius * height_radius_ratio;

  // prepare polygon parameters if needed
  const int max_verts = 8;
  int       n_verts = 6;
  float     thetas[8];
  float     r_verts[8];

  if (shape_type == 1 ||
      shape_type == 3) // SCATTER_SHAPE_POLYGON or SCATTER_SHAPE_PYRAMID
  {
    uint rng = wang_hash(seed ^ (uint)(item_id * 1999 + 17));
    n_verts = 5 + (int)(rand(&rng) * 3.99f);
    if (n_verts > max_verts) n_verts = max_verts;
    float base_rot = rand(&rng) * 2.0f * M_PI_F;

    for (int k = 0; k < n_verts; ++k)
    {
      thetas[k] = base_rot + (2.0f * M_PI_F * (float)k) / (float)n_verts;
      r_verts[k] = radius * (0.65f + 0.35f * rand(&rng));
    }
  }

  for (int iy = iy_min; iy <= iy_max; ++iy)
  {
    float wy = ymin + ((float)iy + 0.5f) * py;
    float dy = wy - cy;

    for (int ix = ix_min; ix <= ix_max; ++ix)
    {
      float wx = xmin + ((float)ix + 0.5f) * px;
      float dx = wx - cx;
      float dist = sqrt(dx * dx + dy * dy);

      float shape_val = 0.0f;

      if (shape_type == 0) // SCATTER_SHAPE_DISK
      {
        if (dist < radius)
        {
          float u = dist / radius;
          shape_val = sqrt(max(0.0f, 1.0f - u * u));
        }
      }
      else if (shape_type == 2) // SCATTER_SHAPE_CONE
      {
        if (dist < radius)
        {
          float u = dist / radius;
          shape_val = 1.0f - u;
        }
      }
      else if (shape_type == 4) // SCATTER_SHAPE_SMOOTH_DOME
      {
        if (dist < radius)
        {
          float u = dist / radius;
          shape_val = 0.5f * (1.0f + cos(u * M_PI_F));
        }
      }
      else if (shape_type == 1 ||
               shape_type == 3) // SCATTER_SHAPE_POLYGON or PYRAMID
      {
        if (dist < 1e-7f)
        {
          shape_val = 1.0f;
        }
        else
        {
          float phi = atan2(dy, dx);
          if (phi < 0.0f) phi += 2.0f * M_PI_F;

          float r_poly = radius;
          for (int k = 0; k < n_verts; ++k)
          {
            int   next_k = (k + 1) % n_verts;
            float t1 = thetas[k];
            float t2 = thetas[next_k];
            while (t1 < 0.0f)
              t1 += 2.0f * M_PI_F;
            while (t1 >= 2.0f * M_PI_F)
              t1 -= 2.0f * M_PI_F;
            while (t2 < 0.0f)
              t2 += 2.0f * M_PI_F;
            while (t2 >= 2.0f * M_PI_F)
              t2 -= 2.0f * M_PI_F;

            bool in_sector = (t1 < t2) ? (phi >= t1 && phi <= t2)
                                       : (phi >= t1 || phi <= t2);
            if (in_sector)
            {
              float d_theta = t2 - t1;
              if (d_theta < 0.0f) d_theta += 2.0f * M_PI_F;
              float d_phi1 = phi - t1;
              if (d_phi1 < 0.0f) d_phi1 += 2.0f * M_PI_F;
              float d_phi2 = t2 - phi;
              if (d_phi2 < 0.0f) d_phi2 += 2.0f * M_PI_F;

              float r1 = r_verts[k];
              float r2 = r_verts[next_k];
              float denom = r2 * sin(d_phi2) + r1 * sin(d_phi1);
              if (denom > 1e-6f)
              {
                r_poly = (r1 * r2 * sin(d_theta)) / denom;
              }
              break;
            }
          }

          if (dist < r_poly)
          {
            float u = dist / r_poly;
            if (shape_type == 1) // POLYGON
            {
              shape_val = sqrt(max(0.0f, 1.0f - u * u));
            }
            else // PYRAMID
            {
              shape_val = 1.0f - u;
            }
          }
        }
      }

      if (shape_val > 0.0f)
      {
        int   idx = linear_index(ix, iy, nx);
        float dz = h_max * shape_val;
        if (dz > 0.0f)
        {
          heightmap[idx] += dz;
        }
        if (has_rock_map != 0)
        {
          rock_map[idx] += shape_val;
        }
      }
    }
  }
}
)""
