#ifndef SITE_H
#define SITE_H

#include "global.hpp"
#include <vector>
#include <algorithm>
#include <iostream>

class Site {
public:
Site(int,int,int,int);
bool step_is_valid(int);
// to change the state of the Site
void set_neighbors(std::vector<Site *>);
void set_neighbors_occ(std::vector<bool*>);
void set_blocking_sites_occ(std::vector<std::vector<bool*> >);
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
Site * move(int);
// debug
void db_print_vector(std::vector<Site *>);
void db_print_properties();
// TODO: check if this needs to be public...

protected:
int m_id;
int m_x;
int m_y;
bool m_is_empty;
int m_type;
std::vector<Site *> m_neighbor_sites;
std::vector<bool *> m_neighbor_is_empty;
std::vector< std::vector <bool*> > m_blocking_sites_occ;
};

#endif
