R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

float schott_get_slope(read_only image2d_t z,
                       const sampler_t     sampler,
                       const int2          p,
                       const int2          q)
{
  float d = length((float2)((p.x - q.x), (p.y - q.y)));
  return (read_imagef(z, sampler, p).x - read_imagef(z, sampler, q).x) / d;
}

float find_steepest_downslope_neighbor(read_only image2d_t z,
                                       const sampler_t     sampler,
                                       const int2          p,
                                       int2               *ptr_q)
{
  float slope_max = 0.f;
  (*ptr_q) = p;

  for (int r = -1; r <= 1; r++)
    for (int s = -1; s <= 1; s++)
      if (!(r == 0 && s == 0))
      {
        int2  pn = (int2)(p.x + r, p.y + s);
        float slope = schott_get_slope(z, sampler, p, pn);

        if (slope > slope_max)
        {
          slope_max = slope;
          (*ptr_q) = pn;
        }
      }

  return slope_max;
}

float schott_get_weight(read_only image2d_t z,
                        const sampler_t     sampler,
                        const int2          p,
                        const int2          q,
                        const float         flow_routing_exponent)
{
  float weight = 0.f;
  float total_weight = 0.f;

  for (int r = -1; r <= 1; r++)
    for (int s = -1; s <= 1; s++)
      if (!(r == 0 && s == 0))
      {
        int2  pn = (int2)(p.x + r, p.y + s);
        float slope = schott_get_slope(z, sampler, p, pn);
        if (slope > 0.f)
        {
          float sp = pow_float(slope, flow_routing_exponent);
          total_weight += sp;

          if (pn.x == q.x && pn.y == q.y) weight = sp;
        }
      }

  if (total_weight == 0.f)
    return 0.f;
  else
    return weight / total_weight;
}

float schott_get_flow(read_only image2d_t z,
                      read_only image2d_t field,
                      const sampler_t     sampler,
                      const int2          p,
                      const float         flow_routing_exponent)
{
  float field_out = 0.f;

  for (int r = -1; r <= 1; r++)
    for (int s = -1; s <= 1; s++)
      if (!(r == 0 && s == 0))
      {
        int2 q = (int2)(p.x + r, p.y + s);
        if (read_imagef(z, sampler, p).x <= read_imagef(z, sampler, q).x)
        {
          float weight = schott_get_weight(z,
                                           sampler,
                                           q,
                                           p,
                                           flow_routing_exponent);
          field_out += read_imagef(field, sampler, q).x * weight;
        }
      }

  return field_out;
}

void kernel hydraulic_schott(read_only image2d_t  z,
                             read_only image2d_t  flow,
                             read_only image2d_t  sediment,
                             read_only image2d_t  talus,
                             write_only image2d_t z_new,
                             write_only image2d_t flow_new,
                             write_only image2d_t sediment_new,
                             const int            nx,
                             const int            ny,
                             const float          c_erosion,
                             const float          c_thermal,
                             const float          c_deposition,
                             const float          flow_acc_exponent,
                             const float          flow_acc_exponent_depo,
                             const float          flow_routing_exponent,
                             const int            erosion_it,
                             const int            thermal_it,
                             const int            it)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  // neighbor
  int2  nbrs;
  float slope_max = find_steepest_downslope_neighbor(z, sampler, g, &nbrs);
  float speed = clamp(slope_max * slope_max, 0.f, 1.f);

  float z_val = read_imagef(z, sampler, g).x;

  if (it % 10 < erosion_it)
  {
    // --- hydraulic erosion

    float spe = c_erosion * min(3.f,
                                pow_float(read_imagef(flow, sampler, g).x,
                                          flow_acc_exponent) *
                                    speed);

    float z_steepest = read_imagef(z, sampler, nbrs).x;

    z_val = max(z_steepest, z_val - spe);
  }
  else if (it % 10 < thermal_it)
  {
    // --- thermal erosion

    float talus_val = read_imagef(talus, sampler, g).x;

    int up = 0;
    int down = 0;

    for (int r = -1; r <= 1; r++)
      for (int s = -1; s <= 1; s++)
        if (!(r == 0 && s == 0))
        {
          int2  q = (int2)(g.x + r, g.y + s);
          float slope = schott_get_slope(z, sampler, g, q);

          if (slope > talus_val)
            down++;
          else if (slope < -talus_val)
            up++;
        }

    z_val += c_thermal * (float)(up - down) * talus_val;
  }
  else
  {
    // --- deposition

    float spe = pow_float(read_imagef(flow, sampler, g).x,
                          flow_acc_exponent_depo) *
                speed;

    float sed = read_imagef(sediment, sampler, g).x;
    sed = schott_get_flow(z, sediment, sampler, g, flow_routing_exponent);

    if (sed > spe)
    {
      float deposit = min(sed, c_deposition * (sed - spe));
      z_val += deposit;
      sed -= deposit;
    }

    sed += c_deposition * spe;

    write_imagef(sediment_new, g, sed);
  }

  write_imagef(z_new, g, z_val);

  float flow_val = 1.f +
                   schott_get_flow(z, flow, sampler, g, flow_routing_exponent);
  write_imagef(flow_new, g, flow_val);
}

void kernel hydraulic_schott_erosion(read_only image2d_t  z,
                                     read_only image2d_t  flow,
                                     read_only image2d_t  moisture_map,
                                     write_only image2d_t z_new,
                                     write_only image2d_t flow_new,
                                     const int            nx,
                                     const int            ny,
                                     const float          c_erosion,
                                     const float          flow_acc_exponent,
                                     const float          flow_routing_exponent)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  // neighbor
  int2  nbrs;
  float slope_max = find_steepest_downslope_neighbor(z, sampler, g, &nbrs);
  float speed = clamp(slope_max * slope_max, 0.f, 1.f);

  float z_val = read_imagef(z, sampler, g).x;

  // --- hydraulic erosion

  float spe = c_erosion * min(3.f,
                              pow_float(read_imagef(flow, sampler, g).x,
                                        flow_acc_exponent) *
                                  speed);

  float z_steepest = read_imagef(z, sampler, nbrs).x;
  z_val = max(z_steepest, z_val - spe);
  write_imagef(z_new, g, z_val);

  float flow_val = read_imagef(moisture_map, sampler, g).x +
                   schott_get_flow(z, flow, sampler, g, flow_routing_exponent);
  write_imagef(flow_new, g, flow_val);
}
)""
