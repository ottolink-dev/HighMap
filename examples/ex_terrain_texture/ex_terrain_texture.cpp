#include <string>
#include <vector>

#include "highmap.hpp"

// Helper function to horizontally stack multiple textures into a single banner
static hmap::Texture hstack_textures(const std::vector<hmap::Texture> &textures)
{
  if (textures.empty())
  {
    return hmap::Texture();
  }

  int total_width = 0;
  int height = textures[0].shape.y;
  int num_channels = textures[0].num_channels();

  for (const auto &tex : textures)
  {
    total_width += tex.shape.x;
  }

  hmap::Texture banner({total_width, height}, num_channels);

  int x_offset = 0;
  for (const auto &tex : textures)
  {
    for (int c = 0; c < num_channels; ++c)
    {
      const auto &src = tex.channels[c];
      auto       &dst = banner.channels[c];

      for (int j = 0; j < tex.shape.y; ++j)
      {
        for (int i = 0; i < tex.shape.x; ++i)
        {
          dst(x_offset + i, j) = src(i, j);
        }
      }
    }

    x_offset += tex.shape.x;
  }

  return banner;
}

// Preset color sets
static hmap::TerrainTextureParams preset_alpine()
{
  hmap::TerrainTextureParams p;

  p.grass = {0.28f, 0.36f, 0.22f};
  p.dry = {0.52f, 0.49f, 0.36f};
  p.soil = {0.45f, 0.38f, 0.30f};
  p.sediment = {0.68f, 0.65f, 0.58f};
  p.rock = {0.42f, 0.41f, 0.40f};
  p.ridge = {0.62f, 0.62f, 0.62f};
  p.snow = {0.95f, 0.96f, 0.98f};

  p.snow_line = 0.70f;
  p.snow_amount = 0.85f;
  p.rock_talus = 0.0035f;
  p.rock_edge = 0.16f;
  p.strata = 0.40f;
  p.streaks = 0.85f;

  return p;
}

static hmap::TerrainTextureParams preset_everest()
{
  hmap::TerrainTextureParams p;

  p.grass = {0.35f, 0.34f, 0.30f};
  p.dry = {0.45f, 0.42f, 0.38f};
  p.soil = {0.36f, 0.33f, 0.29f};
  p.sediment = {0.60f, 0.58f, 0.55f};
  p.rock = {0.28f, 0.27f, 0.28f};
  p.ridge = {0.50f, 0.49f, 0.52f};
  p.snow = {0.96f, 0.97f, 1.00f};

  p.snow_line = 0.45f;
  p.snow_amount = 1.0f;
  p.rock_talus = 0.0025f;
  p.rock_edge = 0.18f;
  p.strata = 0.50f;
  p.streaks = 0.90f;

  return p;
}

static hmap::TerrainTextureParams preset_grassland()
{
  hmap::TerrainTextureParams p;

  p.grass = {0.26f, 0.48f, 0.18f};
  p.dry = {0.58f, 0.60f, 0.28f};
  p.soil = {0.48f, 0.38f, 0.26f};
  p.sediment = {0.66f, 0.58f, 0.42f};
  p.rock = {0.50f, 0.46f, 0.38f};
  p.ridge = {0.65f, 0.60f, 0.50f};
  p.snow = {0.94f, 0.95f, 0.97f};

  p.snow_line = 1.20f;
  p.snow_amount = 0.0f;
  p.rock_talus = 0.0055f;
  p.rock_edge = 0.20f;
  p.strata = 0.20f;
  p.flow = 0.9f;

  return p;
}

static hmap::TerrainTextureParams preset_iceland_black_gray()
{
  hmap::TerrainTextureParams p;

  p.grass = {0.18f, 0.18f, 0.18f};
  p.dry = {0.28f, 0.27f, 0.26f};
  p.soil = {0.14f, 0.13f, 0.13f};
  p.sediment = {0.35f, 0.34f, 0.33f};
  p.rock = {0.10f, 0.10f, 0.11f};
  p.ridge = {0.42f, 0.42f, 0.43f};
  p.snow = {0.92f, 0.94f, 0.97f};

  p.snow_line = 0.80f;
  p.snow_amount = 0.45f;
  p.rock_talus = 0.0030f;
  p.rock_edge = 0.14f;
  p.strata = 0.30f;
  p.cavity = 0.8f;

  return p;
}

static hmap::TerrainTextureParams preset_iceland_black_green()
{
  hmap::TerrainTextureParams p;

  p.grass = {0.22f, 0.46f, 0.15f};
  p.dry = {0.34f, 0.42f, 0.18f};
  p.soil = {0.15f, 0.15f, 0.14f};
  p.sediment = {0.24f, 0.24f, 0.22f};
  p.rock = {0.11f, 0.11f, 0.12f};
  p.ridge = {0.30f, 0.30f, 0.32f};
  p.snow = {0.93f, 0.95f, 0.98f};

  p.snow_line = 0.85f;
  p.snow_amount = 0.35f;
  p.rock_talus = 0.0035f;
  p.rock_edge = 0.15f;
  p.strata = 0.25f;

  return p;
}

static hmap::TerrainTextureParams preset_atoll_island()
{
  hmap::TerrainTextureParams p;

  p.grass = {0.12f, 0.52f, 0.22f};
  p.dry = {0.78f, 0.76f, 0.45f};
  p.soil = {0.84f, 0.79f, 0.65f};
  p.sediment = {0.90f, 0.87f, 0.75f};
  p.rock = {0.62f, 0.58f, 0.52f};
  p.ridge = {0.75f, 0.72f, 0.66f};
  p.snow = {0.98f, 0.98f, 0.98f};

  p.snow_line = 1.50f;
  p.snow_amount = 0.0f;
  p.rock_talus = 0.0050f;
  p.flow = 0.95f;
  p.deposition = 0.70f;

  return p;
}

static hmap::TerrainTextureParams preset_red_canyon()
{
  hmap::TerrainTextureParams p;

  p.grass = {0.46f, 0.38f, 0.25f};
  p.dry = {0.68f, 0.45f, 0.28f};
  p.soil = {0.62f, 0.32f, 0.18f};
  p.sediment = {0.82f, 0.62f, 0.42f};
  p.rock = {0.66f, 0.28f, 0.16f};
  p.ridge = {0.85f, 0.52f, 0.30f};
  p.snow = {0.94f, 0.92f, 0.90f};

  p.snow_line = 1.50f;
  p.snow_amount = 0.0f;
  p.rock_talus = 0.0030f;
  p.rock_edge = 0.18f;
  p.strata = 0.65f;
  p.streaks = 0.75f;

  return p;
}

static hmap::TerrainTextureParams preset_desert_dunes()
{
  hmap::TerrainTextureParams p;

  p.grass = {0.74f, 0.62f, 0.42f};
  p.dry = {0.82f, 0.68f, 0.46f};
  p.soil = {0.70f, 0.56f, 0.38f};
  p.sediment = {0.88f, 0.76f, 0.52f};
  p.rock = {0.65f, 0.52f, 0.36f};
  p.ridge = {0.92f, 0.82f, 0.60f};
  p.snow = {0.96f, 0.94f, 0.92f};

  p.snow_line = 1.50f;
  p.snow_amount = 0.0f;
  p.rock_talus = 0.0060f;
  p.strata = 0.25f;
  p.streaks = 0.60f;

  return p;
}

static hmap::TerrainTextureParams preset_alien_sulfur()
{
  hmap::TerrainTextureParams p;

  p.grass = {0.72f, 0.68f, 0.14f};
  p.dry = {0.58f, 0.52f, 0.12f};
  p.soil = {0.38f, 0.22f, 0.15f};
  p.sediment = {0.85f, 0.80f, 0.20f};
  p.rock = {0.26f, 0.20f, 0.24f};
  p.ridge = {0.88f, 0.62f, 0.15f};
  p.snow = {0.92f, 0.90f, 0.75f};

  p.snow_line = 1.50f;
  p.snow_amount = 0.0f;
  p.rock_talus = 0.0035f;
  p.rock_edge = 0.16f;
  p.strata = 0.45f;
  p.streaks = 0.85f;

  return p;
}

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2    shape = {1024, 1024};
  std::uint32_t seed = 0;

  // Generate a procedural mountainous terrain
  hmap::Array z = hmap::gpu::mountain_stump(shape, seed);

  // Apply some hydraulic erosion to generate deposition and flow channels
  hmap::Array deposition(shape, 0.f);
  hmap::gpu::hydraulic_particle(z, 200000, seed, nullptr, nullptr, &deposition);

  // Build the list of color presets
  struct PresetEntry
  {
    std::string                name;
    hmap::TerrainTextureParams params;
  };

  std::vector<PresetEntry> presets = {
      {"alpine", preset_alpine()},
      {"everest", preset_everest()},
      {"grassland", preset_grassland()},
      {"iceland_black_gray", preset_iceland_black_gray()},
      {"iceland_black_green", preset_iceland_black_green()},
      {"atoll_island", preset_atoll_island()},
      {"red_canyon", preset_red_canyon()},
      {"desert_dunes", preset_desert_dunes()},
      {"alien_sulfur", preset_alien_sulfur()},
  };

  std::vector<hmap::Texture> textures;
  textures.reserve(presets.size());

  for (const auto &entry : presets)
  {
    hmap::Texture tex = hmap::terrain_texture(z,
                                              &deposition,
                                              nullptr,
                                              entry.params);
    tex.to_png("ex_terrain_texture_" + entry.name + ".png");
    textures.push_back(std::move(tex));
  }

  // Horizontal banner combining all presets
  hmap::Texture banner = hstack_textures(textures);
  banner.to_png("ex_terrain_texture_banner.png");

  z.dump();

  return 0;
}
