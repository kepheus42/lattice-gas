#ifndef SITE_H
#define SITE_H

#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>

class Site {
public:
Site(int,int,int,int);
bool step_is_invalid(int);
bool get_blocking_site_state(int,int);
std::vector<unsigned int> blocking_site_corr();
std::vector<std::vector<Site *> > get_blocking_sites();
std::vector<int> get_site_correlation();
// to change the state of the Site
void set_neighbor_sites(std::vector<Site *>);
void set_blocking_sites(std::vector<std::vector<Site *> >);
void swap_state();
void set_not_empty();
// to get the state of the Site
bool is_empty();
bool * get_state_ptr();
int get_id();
int get_x();
int get_y();
Site * get_neighbor_by_dir(int);
// move from one site to a neighbor, specified by dir
Site * jump_in_direction(int);
// debug
void db_print_vector(std::vector<Site *>);
void db_print_properties();

void print_neighbor_sites();
void print_blocking_sites();
// TODO: check if this needs to be public...

protected:
int m_id;
int m_x;
int m_y;
bool m_is_empty;
int m_type;
int m_blocking_sites_per_dir;
std::vector<Site *> m_neighbor_sites;
std::vector< std::vector <Site *> > m_blocking_sites;
std::vector<int> m_powers_of_two;
};

#endif
