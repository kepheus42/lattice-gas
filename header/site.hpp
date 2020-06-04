#ifndef SITE_H
#define SITE_H

#include <vector>
#include <algorithm>

class Site {
public:
Site(int,int,int);
void set_neighbors_1x1(std::vector<Site *>);
void set_neighbors_2x2(std::vector<Site *>);
// to change the state of the Site
void set_empty();
void set_occupied();
// to get the state of the Site
bool is_empty();
bool is_blocked();
//
virtual Site * get_neighbor_by_dir(int);
//
virtual std::vector<int> get_correlations();
//
int get_x();
int get_y();

private:
int m_id;
int m_x;
int m_y;
bool m_is_empty;
std::vector<Site *> m_neighbors_1x1;
std::vector<Site *> m_neighbors_2x2;
};

class Site_1x1 : public Site {
public:
Site_1x1(int,int,int);
Site * get_neighbor_by_dir(int);
bool step_is_valid(int);
std::vector<int> get_correlations();
private:
};

class Site_2x2 : public Site {
public:
Site_2x2(int,int,int);
Site * get_neighbor_by_dir(int);
bool step_is_valid(int);
std::vector<int> get_correlations();
private:
};


#endif
