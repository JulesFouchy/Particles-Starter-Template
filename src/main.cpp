#include <iostream>
#include "glm/ext/quaternion_common.hpp"
#include "glm/ext/scalar_constants.hpp"
#include "glm/geometric.hpp"
#include "opengl-framework/opengl-framework.hpp"
#include "utils.hpp"

float easeInOut(float x, float power)
{
    if (x < 0.5)
    {
        return 0.5 * pow(2 * x, power);
    }
    else
    {
        return 1 - 0.5 * pow(2 * (1 - x), power);
    }
}

struct Wall {
    glm::vec2 start{};
    glm::vec2 end{};
};

struct Segment {
    glm::vec2 start{};
    glm::vec2 end{};
};

bool nearly_equal(float a, float b, float epsilon = 1e-6f)
{
    return std::fabs(a - b) < epsilon;
}

bool is_invertible(const glm::mat2& mat, float epsilon = 1e-6f)
{
    return !nearly_equal(glm::determinant(mat), 0.0f, epsilon);
}

struct Circle {
    glm::vec2 center;
    float     radius;
};

std::optional<glm::vec2> intersection(Segment const& segment, Circle const& circle)
{
    glm::vec2 const o      = segment.start;
    glm::vec2 const d      = segment.end - segment.start;
    glm::vec2 const center = circle.center;
    float const     r      = circle.radius;

    float const a = glm::dot(d, d);
    float const b = 2 * glm::dot(o - center, d);
    float const c = glm::dot(o - center, o - center) - r * r;

    float const delta = b * b - 4 * a * c;
    if (delta <= 0)
        return std::nullopt;

    float const t1 = (-b - std::sqrt(delta)) / 2.f / a;
    float const t2 = (-b + std::sqrt(delta)) / 2.f / a;

    if (0 <= t1 && t1 <= 1 && 0 <= t2 && t2 <= 1)
        return o + std::min(t1, t2) * d;

    if (0 <= t1 && t1 <= 1)
        return o + t1 * d;

    if (0 <= t2 && t2 <= 1)
        return o + t2 * d;

    return std::nullopt;
}

std::optional<glm::vec2> intersection(Segment const& s1, Segment const& s2)
{
    glm::vec2 const o1 = s1.start;
    glm::vec2 const o2 = s2.start;
    glm::vec2 const d1 = s1.end - s1.start;
    glm::vec2 const d2 = s2.end - s2.start;

    glm::mat2x2 const mat = glm::mat2x2{d1, -d2};
    // if (!is_invertible(mat))
    //     return std::nullopt;

    glm::vec2 const solutions = glm::inverse(mat) * (o2 - o1);
    if (0 <= solutions.x && solutions.x <= 1 && 0 <= solutions.y && solutions.y <= 1)
        return o1 + solutions.x * d1;

    return std::nullopt;
}

struct Particle {
    glm::vec2 position{
        utils::rand(-gl::window_aspect_ratio(), +gl::window_aspect_ratio()),
        utils::rand(-1.f, +1.f),
    };

    glm::vec2 velocity;

    float mass{utils::rand(1.f, 2.f)};

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
        return glm::mix(start_color, end_color, easeInOut(relative_age(), 4.f));
    }

    float radius() const
    {
        return std::min(lifespan - age, 2.f) / 2.f * 0.03f;
    }

    float relative_age() const
    {
        return age / lifespan;
    }

    Particle()
    {
        float const initial_angle = utils::rand(0.f, 2.f * glm::pi<float>());

        velocity = {
            0.2f * std::cos(initial_angle),
            0.2f * std::sin(initial_angle),
        };
    }
};

int main()
{
    gl::init("Particules!");
    gl::maximize_window();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    std::vector<Particle> particles(100);

    std::vector<Wall> walls;
    walls.push_back(Wall{{-0.3f, -0.4f}, {0.1f, 0.8f}});
    walls.push_back(Wall{{0.f, 0.f}, {0.3f, 0.3f}});

    std::vector<Circle> circles;
    circles.push_back(Circle{.center = {0.3f, 0.f}, .radius = 0.2f});
    circles.push_back(Circle{.center = {0.7f, -0.5f}, .radius = 0.4f});
    circles.push_back(Circle{.center = {-0.4f, 0.f}, .radius = 0.1f});

    for (auto const& circle : circles)
    {
        std::erase_if(particles, [&](Particle const& particle) { return glm::distance(particle.position, circle.center) < circle.radius; });
    }

    std::optional<glm::vec2> place_wall{};

    gl::set_events_callbacks({gl::EventsCallbacks{
        .on_mouse_pressed = [&](gl::MousePressedEvent const& e) {
            auto const w = gl::window_width_in_screen_coordinates();
            auto const h = gl::window_height_in_screen_coordinates();

            auto pos = (e.position - glm::vec2{w, h} * 0.5f) / (float)gl::window_height_in_screen_coordinates() * 2.f;
            pos.y *= -1.f;
            if (place_wall.has_value())
            {
                walls.push_back(Wall{*place_wall, pos});
                place_wall.reset();
                std::cout << walls.size() << '\n';
            }
            else
            {
                place_wall = pos;
            }
        },
    }});

    while (gl::window_is_open())
    {
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);

        // for (auto const& circle : circles)
        // {
        //     // auto const segment1 = Segment{{-1.f, 0.f}, {1.f, 0.f}};
        //     auto const segment2 = Segment{{0.f, -1.f}, gl::mouse_position()};
        //     // utils::draw_line(segment1.start, segment1.end, 0.01f, glm::vec4{1.f});
        //     utils::draw_line(segment2.start, segment2.end, 0.01f, glm::vec4{1.f});

        //     auto const inter = intersection(segment2, circle);
        //     if (inter)
        //         utils::draw_disk(*inter, 0.05f, glm::vec4{0., 1.f, 0.f, 1.f});
        // }

        for (auto& particle : particles)
        {
            particle.age += gl::delta_time_in_seconds();

            auto forces = glm::vec2{0.f};

            // Gravity
            // forces += glm::vec2{0.f, -1.f} * particle.mass;

            // Air friction
            // forces += -particle.velocity * 1.f;

            // Follow mouse
            // forces += (gl::mouse_position() - particle.position);

            particle.velocity += forces / particle.mass * gl::delta_time_in_seconds();

            auto const old_pos = particle.position;
            particle.position += particle.velocity * gl::delta_time_in_seconds();
            auto const new_pos = particle.position;

            auto tmpwalls = walls;
            tmpwalls.push_back(Wall{{-gl::window_aspect_ratio(), -1.f}, {+gl::window_aspect_ratio(), -1.f}});
            tmpwalls.push_back(Wall{{+gl::window_aspect_ratio(), -1.f}, {+gl::window_aspect_ratio(), +1.f}});
            tmpwalls.push_back(Wall{{+gl::window_aspect_ratio(), +1.f}, {-gl::window_aspect_ratio(), +1.f}});
            tmpwalls.push_back(Wall{{-gl::window_aspect_ratio(), +1.f}, {-gl::window_aspect_ratio(), -1.f}});

            for (auto const& wall : tmpwalls)
            {
                auto const inter = intersection({old_pos, new_pos}, {wall.end, wall.start});
                if (!inter)
                    continue;
                auto const wall_dir    = glm::normalize(wall.end - wall.start);
                auto const wall_normal = glm::vec2{-wall_dir.y, wall_dir.x};
                particle.velocity      = glm::reflect(particle.velocity, wall_normal);

                auto desired_dist   = glm::distance(old_pos, new_pos);
                auto dist_to_wall   = glm::distance(old_pos, *inter);
                auto remaining_dist = desired_dist - dist_to_wall;
                particle.position   = *inter + glm::normalize(particle.velocity) * remaining_dist;
            }
            for (auto const& circle : circles)
            {
                auto const inter = intersection(Segment{old_pos, new_pos}, circle);
                if (!inter)
                    continue;
                auto const normal = glm::normalize(*inter - circle.center);
                particle.velocity = glm::reflect(particle.velocity, normal);

                auto desired_dist   = glm::distance(old_pos, new_pos);
                auto dist_to_wall   = glm::distance(old_pos, *inter);
                auto remaining_dist = desired_dist - dist_to_wall;
                particle.position   = *inter + glm::normalize(particle.velocity) * remaining_dist;
            }
        }

        std::erase_if(particles, [&](Particle const& particle) { return particle.age > particle.lifespan; });

        for (auto const& particle : particles)
            utils::draw_disk(particle.position, particle.radius(), glm::vec4{particle.color(), 1.f});

        for (auto const& wall : walls)
            utils::draw_line(wall.start, wall.end, 0.01f, glm::vec4{1.f});

        for (auto const& circle : circles)
            utils::draw_disk(circle.center, circle.radius, glm::vec4{1.f, 1., 1., 0.3f});
    }
}