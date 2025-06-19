#include "glm/ext/scalar_constants.hpp"
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

struct Particle {
    glm::vec2 position{
        utils::rand(-gl::window_aspect_ratio(), +gl::window_aspect_ratio()),
        utils::rand(-1.f, +1.f),
    };

    glm::vec2 velocity;

    float mass{utils::rand(1.f, 2.f)};

    float age{0.f};
    float lifespan{utils::rand(3.f, 5.f)};

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

    while (gl::window_is_open())
    {
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);

        for (auto& particle : particles)
        {
            particle.age += gl::delta_time_in_seconds();

            auto forces = glm::vec2{0.f};

            // Gravity
            // forces += glm::vec2{0.f, -1.f} * particle.mass;

            // Air friction
            forces += -particle.velocity * 1.f;

            // Follow mouse
            // forces += (gl::mouse_position() - particle.position);

            particle.velocity += forces / particle.mass * gl::delta_time_in_seconds();
            particle.position += particle.velocity * gl::delta_time_in_seconds();
        }

        std::erase_if(particles, [&](Particle const& particle) { return particle.age > particle.lifespan; });

        for (auto const& particle : particles)
            utils::draw_disk(particle.position, particle.radius(), glm::vec4{particle.color(), 1.f});
    }
}