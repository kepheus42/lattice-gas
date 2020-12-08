#ifndef SITE_H
#define SITE_H


#include <vector>
#include <algorithm>
#include <iostream>

class Site {
public:
Site(int,int,int,int);
bool step_is_valid(int);
bool blocking_site_state(int,int);
std::vector<unsigned int> blocking_site_corr();
std::vector<std::vector<bool *> > get_blocking_sites();
// to change the state of the Site
void set_neighbor_sites(std::vector<Site *>);
void set_blocking_sites(std::vector<std::vector<bool *> >);
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
std::vector< std::vector <bool *> > m_blocking_sites;
};

#endif
