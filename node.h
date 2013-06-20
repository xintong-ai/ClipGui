#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

#define X 400
#define Y 400

#define BG  GetNamedColor("black")
#define WORK  GetNamedColor("white")
#define S  GetNamedColor("cyan")
#define C  GetNamedColor("green")
#define POLY  GetNamedColor("red")
#include <stdbool.h>

enum SIDE{UNCLEAR,IN,OUT,ON};

typedef struct _node
{
  int x, y;
  struct _node *next;
  struct _node *prev;
  struct _node *nextPoly;   /* pointer to the next polygon */
  struct _node *neighbor;   /* the coresponding intersection point */
  int intersect;            /* 1 if an intersection point, 0 otherwise */
  int entry;                /* 1 if an entry point, 0 otherwise */
  int visited;              /* 1 if the node has been visited, 0 otherwise */
  float alpha;              /* intersection point placemet */
  SIDE side;
  _node()
  {
      entry = 1;
      side = UNCLEAR;
  }
} node;
