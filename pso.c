#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

typedef struct {
    double x;
    double y;
    double rotation;
} coordinate_t;

typedef struct {
    coordinate_t lower;
    coordinate_t upper;
} space_t;

typedef struct {
    space_t search_space;
    coordinate_t actual_coord;
} context_t;

typedef struct {
    coordinate_t position;
    coordinate_t velocity;

    coordinate_t best_position;
    double best_position_cost;
} particle_t;

typedef struct {
    int swarm_size;
    double omega;
    double phi_p;
    double phi_g;
} swarm_params_t;

int run_pso(context_t *context, swarm_params_t swarm_params, coordinate_t *best_position, double *best_position_cost);
double calculate_cost(coordinate_t *coord, coordinate_t *actual);

int main(int argc, char **argv) {
    srand(time(NULL));

    space_t search_space = {
        .lower = { 0, 0, 0 },
        .upper = { 999999, 999999, M_PI * 2 }
    };

    context_t context = {
        .search_space = search_space,
        .actual_coord = { 2500, 400, M_PI }
    };

    swarm_params_t swarm_params = {
        .swarm_size = 100,
        .omega = 0.5,
        .phi_p = 0.5,
        .phi_g = 0.5
    };

    coordinate_t best_position;
    double best_position_cost;

    int count = run_pso(&context, swarm_params, &best_position, &best_position_cost);

    fprintf(stdout, "Iterations: %d, best pos: (%.2f, %.2f, %.2f)\n", count, best_position.x, best_position.y, best_position.rotation);
}

/**
 * Gets a random value between the lower and upper bounds
 */
int random_value(const int lower, const int upper, bool include_negatives) {
    int random_val = rand() % (upper - lower);

    if (include_negatives) {
        return random_val * 2 - (upper - lower) + 1;
    }

    return random_val + lower;
}

/**
 * Gets a dimension's reference based on the specified index
 */
double *get_coord(coordinate_t *coord, const int index) {
    switch (index) {
        case 0: return &coord->x;
        case 1: return &coord->y;
        case 2: return &coord->rotation;
    }
    return NULL;
}

/**
 * Generates a random coordinate with each dimension initialized to a random value within the specified bounds
 */
coordinate_t random_coord(space_t *bounds, bool include_negatives) {
    coordinate_t coord;

    // Assuming each dimension is a double, we can loop through them by figuring out how many dimensions there are
    for (int c = 0; c < sizeof(coordinate_t) / sizeof(double); c++) {
        *get_coord(&coord, c) = random_value(
            *get_coord(&bounds->lower, c),
            *get_coord(&bounds->upper, c),
            include_negatives
        );
    }

    return coord;
}

/**
 * Sets the particles to random positions and velocities which are contained inside the search space
 */
void fill_particles(context_t *context, particle_t particles[], const int size, coordinate_t *best_position, double *best_position_cost) {
    *best_position_cost = INFINITY;

    for (int i = 0; i < size; i++) {
        coordinate_t position = random_coord(&context->search_space, false);
        double cost = calculate_cost(&position, &context->actual_coord);

        particles[i].position = position;
        particles[i].velocity = random_coord(&context->search_space, true);
        particles[i].best_position = position;
        particles[i].best_position_cost = cost;

        if (cost < *best_position_cost) {
            *best_position = position;
            *best_position_cost = cost;
        }
    }
}

/**
 * Calculates the new velocity for a single dimension
 */
double calculate_velocity(swarm_params_t *swarm_params, double velocity, double position, double particle_best_position, double best_position) {
    double random_p = random_value(0, 1000, false) / 1000.0;
    double random_g = random_value(0, 1000, false) / 1000.0;

    return swarm_params->omega * velocity +
        swarm_params->phi_p * random_p * (particle_best_position - position) +
        swarm_params->phi_g * random_g * (best_position - position);
}

/**
 * Updates a particle's velocity and position for all dimensions
 */
void update_particle(particle_t *particle, coordinate_t *best_position, swarm_params_t *swarm_params) {
    // Assuming each dimension is a double, we can loop through them by figuring out how many dimensions there are
    for (int c = 0; c < sizeof(coordinate_t) / sizeof(double); c++) {
        double *particle_velocity = get_coord(&particle->velocity, c);
        double *particle_position = get_coord(&particle->position, c);
        double *particle_best_position = get_coord(&particle->best_position, c);
        double *best_position_val = get_coord(best_position, c);

        *particle_velocity = calculate_velocity(swarm_params,
            *particle_velocity, *particle_position, *particle_best_position, *best_position_val
        );

        // TODO: Should this really be working better?
        // Do we need to keep the velocity and position within 2PI range?
        if (c == 2) {
            *particle_velocity = *particle_velocity / M_PI;
        }

        *particle_position += *particle_velocity;
    }
}

/**
 * Runs particle swarm optimisation with the given swarm parameters
 */
int run_pso(context_t *context, swarm_params_t swarm_params, coordinate_t *best_position, double *best_position_cost) {
    int iterations = 0;

    particle_t particles[swarm_params.swarm_size];
    fill_particles(context, particles, swarm_params.swarm_size, best_position, best_position_cost);

    for (int i = 0; i < 1000; i++) {
        for (int p = 0; p < swarm_params.swarm_size; p++) {
            particle_t *particle = &particles[p];
            coordinate_t *position = &particle->position;

            update_particle(particle, best_position, &swarm_params);
            double cost = calculate_cost(position, &context->actual_coord);

            if (cost < particle->best_position_cost) {
                particle->best_position = *position;
                particle->best_position_cost = cost;

                if (cost < *best_position_cost) {
                    *best_position = *position;
                    *best_position_cost = cost;
                }
            }
        }

        iterations++;
    }

    return iterations;
}

/**
 * Calculates the cost from the specified coordinate to the actual position
 */
double calculate_cost(coordinate_t *coord, coordinate_t *actual) {
    return sqrt(
        pow(fabs(actual->x - coord->x) / 1000.0, 2) +
        pow(fabs(actual->y - coord->y) / 1000.0, 2) +
        pow(fabs(actual->rotation - coord->rotation) / M_PI / 2.0, 2)
    );
}
