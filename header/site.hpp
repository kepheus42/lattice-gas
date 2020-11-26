#ifndef SITE_H
#define SITE_H

#include <vector>
#include <algorithm>

/*
   TODO: rewrite to store all relevant connected sites in one vector:
   0 - (z-1) : nearest neighbors (to jump to)
   z - X : connected sites potentially blocked by species 1 tracers
   X - X : connected sites potentially blocked by species 2 tracers
   . . .
   TODO: store relationship of site to neighboring sites in appropriate structure
   such that next sites of jumps in directions 1-z, and sites to check for occupation in those directions are accessible
   Then func: Site * get_nearest_neighbor_by_dir(int dir) should return the appropriate site to perform a tracer move
   and  func: int step_is_invalid(int dir) should return the correct index for step failure counting, or 0 if step is valid
   TODO: write two maps: dir -> int, and: dir -> vector<int>

 */

class Site {
public:
Site(int,int,int,int);
bool step_is_valid(int);
// to change the state of the Site
void set_neighbors(std::vector<Site *>);
void set_neighbor_occupancy(std::vector<std::vector<bool*> >);
void swap_state();
// to get the state of the Site
bool is_empty();
bool * get_state_ptr();
int get_id();
int get_x();
int get_y();
Site * get_neighbor_by_dir(int);
// move from one site to a neighbor, specified by dir
Site * move_to(int);
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
std::vector<Site *> m_neighbors;
std::vector< std::vector <bool*> > m_neighbor_occupancy;
};

class Site_1x1 : public Site {
public:
Site_1x1(int,int,int);
bool step_is_valid(int);
private:
};

class Site_2x2 : public Site {
public:
Site_2x2(int,int,int);
bool step_is_valid(int);
private:
};


#endif
