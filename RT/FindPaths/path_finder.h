#ifndef __PATH_FINDER_H__
#define __PATH_FINDER_H__

enum PTYPE {GPU, CPU};
    
struct pathfinder{
    PTYPE  type;
    pathfinder(enum PTYPE type_):type(type_){}
};

typedef struct pathfinder *pathfinder_p;

#endif
