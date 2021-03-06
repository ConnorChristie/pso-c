# Particle Swarm Optimisation
This is an implementation of the [PSO algorithm](https://en.wikipedia.org/wiki/Particle_swarm_optimization) written in C.

## Running
Call the following method with your search bounds and swarm params:
```C
coordinate_t run_pso(context_t *context, swarm_params_t swarm_params);
```

## Example Code
```C
space_t search_space = {
    .lower = { 0, 0, 0 },
    .upper = { 10000, 10000, M_PI * 2 }
};

context_t context = {
    .search_space = search_space,
    .actual_coord = { 2500, 400, M_PI }
};

swarm_params_t swarm_params = {
    .swarm_size = 100,
    .iterations = 200,
    .omega = 0.5,
    .phi_p = 0.5,
    .phi_g = 0.5
};

coordinate_t best_position = run_pso(&context, swarm_params);
```