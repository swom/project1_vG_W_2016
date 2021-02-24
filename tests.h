#ifndef TESTS_H
#define TESTS_H

#include<cassert>

#include "simulation.h"

void test_all();

void test_demographic_cycle() noexcept;

void test_demographic_sim();

void test_env_changer();

void test_env_grid_cell();

void test_env_param() noexcept;

void test_environment();

void test_funder_data() noexcept;

void test_funders() noexcept;

void test_funders_success() noexcept;

void test_GRN();

void test_ind_param() noexcept;

void test_individual();

void test_meta_param() noexcept;

void test_phenotype();

void test_pop_param() noexcept;

void test_population() noexcept;

void test_sim_param() noexcept;

void test_simulation();

void test_utilities() noexcept;

#endif // TESTS_H
