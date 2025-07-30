#include <cstddef>
#include <iostream>
#include "exe_path/exe_path.h"
#include "glm/ext/quaternion_common.hpp"
#include "glm/ext/scalar_constants.hpp"
#include "glm/ext/vector_float2.hpp"
#include "glm/geometric.hpp"
#include "opengl-framework/opengl-framework.hpp"
#include "utils.hpp"

//
#include <math.h>
#include <stdio.h>
#include <string.h>

struct Particle {
    glm::vec2 position{};

    glm::vec2 velocity{0.f};

    float age{0.f};
    float lifespan{utils::rand(30.f, 50.f)};

    glm::vec3 start_color{
        utils::rand(0.f, 1.f),
        utils::rand(0.f, 1.f),
        utils::rand(0.f, 1.f),
    };
    glm::vec3 end_color{
        utils::rand(0.f, 1.f),
        utils::rand(0.f, 1.f),
        utils::rand(0.f, 1.f),
    };

    glm::vec3 color() const
    {
        return glm::vec3{0.7f, 0.3f, 0.3f}; // glm::mix(start_color, end_color, easeInOut(relative_age(), 4.f));
    }

    float radius() const
    {
        return 0.02f / 2.f; // std::min(lifespan - age, 2.f) / 2.f * 0.03f;
    }

    float relative_age() const
    {
        return age / lifespan;
    }
};

glm::vec2 bezier1(glm::vec2 p1, glm::vec2 p2, float t)
{
    return glm::mix(p1, p2, t);
}

glm::vec2 bezier2(glm::vec2 p1, glm::vec2 a, glm::vec2 p2, float t)
{
    return bezier1(bezier1(p1, a, t), bezier1(a, p2, t), t);
}

glm::vec2 bezier3(glm::vec2 p1, glm::vec2 a1, glm::vec2 a2, glm::vec2 p2, float t)
{
    return bezier1(bezier2(p1, a1, a2, t), bezier2(a1, a2, p2, t), t);
}

void draw_parametric(std::function<glm::vec2(float)> const& parametric)
{
    glm::vec2    prevPoint = parametric(0.);
    size_t const nb_points{100};
    for (size_t i = 1; i <= nb_points; ++i)
    {
        glm::vec2 const point = parametric((float)i / (float)nb_points);
        utils::draw_line(prevPoint, point, 0.03, glm::vec4{1.f, 1.f, 1.f, 1.f});
        prevPoint = point;
    }
}

// glm::vec2 closest(std::function<glm::vec2(float)> const& parametric, glm::vec2 const& point)
// {
//     auto const distance = [&](float t) {
//         return glm::distance(parametric(t), point);
//     };
//     float best_t{0.5f};
//     for (int _ = 0; _ < 100; ++_)
//     {
//         float t = 0.5f; // utils::rand(0.f, 1.f);

//         for (int i = 0; i < 500; ++i)
//         {
//             float const deriv = (distance(t + 0.001) - distance(t)) / 0.001;
//             t -= deriv * 0.001f;
//             t = std::clamp(t, 0.f, 1.f);
//         }
//         if (distance(t) < distance(best_t))
//             best_t = t;
//     }
//     return parametric(best_t);
// }
glm::vec2 closest(std::function<glm::vec2(float)> const& parametric, glm::vec2 const& point, float* out_t = nullptr)
{
    auto const distance = [&](float t) {
        return glm::distance(parametric(t), point);
    };
    float t{0.};
    for (size_t i = 1; i < 10; ++i)
    {
        float const tt = (float)i / 9.f;
        if (distance(tt) < distance(t))
            t = tt;
    }

    for (int i = 0; i < 500; ++i)
    {
        float const deriv = (distance(t + 0.001) - distance(t)) / 0.001;
        t -= deriv * 0.001f;
        t = std::clamp(t, 0.f, 1.f);
    }
    if (out_t)
        *out_t = t;
    return parametric(t);
}

int main()
{
    gl::init("Particules!");
    gl::maximize_window();
    // glEnable(GL_BLEND);
    // glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    auto const parametric = [&](float t) {
        // Bezier
        // return bezier3({-.3f, -.3f}, {-0.2f, 0.5f}, {0.6f, -0.5f}, {.8f, .5f}, t);

        // Heart
        t       = 2.f * t * glm::pi<float>();
        float s = sin(t);
        return glm::vec2{
                   16 * s * s * s,
                   13 * cos(t) - 5 * cos(2 * t) - 2 * cos(3 * t) - cos(4 * t)
               }
               / 20.f;
    };

    std::vector<Particle> particles(100);
    // for (auto& particle : particles)
    //     particle.position = parametric(utils::rand(0.f, 1.f));
    for (size_t i = 0; i < particles.size(); ++i)
    {
        float const t         = (float)i / (float)(particles.size() - 1);
        particles[i].position = parametric(t);
        auto const tangent    = parametric(t + 0.001) - parametric(t - 0.001);
        auto const normal     = glm::normalize(glm::vec2{-tangent.y, tangent.x});
        particles[i].velocity = normal * 0.2f;
    }

    for (auto& particle : particles)
    {
        particle.position.x = utils::rand(-1.f, 1.f);
        particle.position.y = 1.f;
        particle.velocity   = glm::vec2{0.f};
    }

    bool start = true; // false;

    gl::set_events_callbacks({gl::EventsCallbacks{
        .on_scroll = [&](auto&&) {
            start = true;
        },
    }});

    while (gl::window_is_open())
    {
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);

        // utils::draw_disk({-.3f, -.3f}, 0.02, {0.5f, 0.f, 0.5f, 1.f});
        // utils::draw_disk({-0.2f, 0.5f}, 0.02, {0.5f, 0.f, 0.5f, 1.f});
        // utils::draw_disk({.8f, .5f}, 0.02, {0.5f, 0.f, 0.5f, 1.f});
        // utils::draw_disk({0.6f, -0.5f}, 0.02, {0.5f, 0.f, 0.5f, 1.f});
        // utils::draw_disk(gl::mouse_position(), 0.02, {0.5f, 0.f, 0.5f, 1.f});
        // draw_parametric([](float t) {
        //     t       = 2.f * t * glm::pi<float>();
        //     float s = sin(t);
        //     return glm::vec2{
        //                16 * s * s * s,
        //                13 * cos(t) - 5 * cos(2 * t) - 2 * cos(3 * t) - cos(4 * t)
        //            }
        //            / 20.f;
        // });
        // draw_parametric([](float t) {
        //     return bezier1({-.3f, -.3f}, gl::mouse_position(), t);
        // });
        // draw_parametric([](float t) {
        //     return bezier2({-.3f, -.3f}, gl::mouse_position(), {.8f, .5f}, t);
        // });
        draw_parametric(parametric);

        // utils::draw_disk(closest(parametric, gl::mouse_position()), 0.02, {0.5f, 0.f, 0.5f, 1.f});
        // utils::draw_disk(gl::mouse_position(), 0.02, {0.5f, 1.f, .5f, 1.f});

        if (start)
        {
            for (auto& particle : particles)
            {
                particle.age += gl::delta_time_in_seconds();

                auto forces = glm::vec2{0.f};

                // Gravity
                forces += glm::vec2{0.f, -.3f};

                // Bezier
                float      t;
                auto const point   = closest(parametric, particle.position, &t);
                auto const tangent = glm::normalize(parametric(t + 0.001) - parametric(t - 0.001));
                auto const normal  = glm::vec2{-tangent.y, tangent.x};
                auto const dist    = std::max(glm::distance(point, particle.position) * 10.f, 0.5f);
                forces += normal / dist / dist;

                // Air friction
                // forces += -particle.velocity * 1.f;

                // Follow mouse
                // forces += (gl::mouse_position() - particle.position);

                particle.velocity += forces * gl::delta_time_in_seconds();
                particle.position += particle.velocity * gl::delta_time_in_seconds();
            }
        }

        // std::erase_if(particles, [&](Particle const& particle) { return particle.age > particle.lifespan; });

        for (auto const& particle : particles)
            utils::draw_disk(particle.position, particle.radius(), glm::vec4{particle.color(), 1.f});

        for (auto& particle : particles)
        {
            if (particle.position.y < -1.3f)
            {
                particle.position.x = utils::rand(-1.f, 1.f);
                particle.position.y = 1.f;
                particle.velocity   = glm::vec2{0.f};
            }
        }
    }
}