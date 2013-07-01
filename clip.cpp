#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include "libsx.h"
#include "node.h"
#include "math.h"
#include "stdlib.h"
#include "float.h"
#include "algorithm"
#include "vector"
#define TRADITIONAL 0
#define EPS 0.00001

node *s=0, *c=0, *root=0;
int DRAW=1, CLIP=1, pS=1, pC=1;
Widget W[8];

struct pt
{
    float x;
    float y;
 //   bool in;    //either inside or intersection point
    float loc;
    pt(float _x, float _y)
    {
        x = _x;
        y = _y;
        loc = -1;
    }
    pt()
    {
        loc = -1;
    }
};

struct trgl
{
    //the first 3 points are the vertex
    //others are reserved forintersection points
    pt p[9];
};

void view_node(node *p)
{
  if(p) printf("%c%c%c (%3d,%3d)  %f    c:%10p n:%10p P:%10p\n",
        p->intersect ? 'I' : ' ',
        p->entry ? 'E' : ' ',
        p->visited ? 'X' : ' ',
        p->x, p->y, p->alpha, p, p->neighbor, p->nextPoly);
  else  puts("NULL");
}

void view(node *p)
{
  node *aux=p;
  puts("");

  if(aux) do
  {
        view_node(aux);
        aux=aux->next;
  }
  while(aux && aux != p);
}

void plot(node *p)
{
  node *aux=p;
  SetColor( WORK);

  if(aux) do
  {
        DrawLine(aux->x-2, aux->y, aux->x+2, aux->y);
        DrawLine(aux->x, aux->y-2, aux->x, aux->y+2);
        aux=aux->next;
  }
  while(aux && aux != p);
}

void deleteNode(node *p)
{
  node *aux, *hold;

  if(hold=p) do
  {
        aux=p;
        p=p->next;
        free(aux);
  }
  while(p && p!=hold);
}

void removeNode(node *p)
{
    p->prev->next = p->next;
    p->next->prev = p->prev;
    free(p);
}

void insert(node *ins, node *first, node *last)
{
  node *aux=first;
  while(aux != last && aux->alpha <= ins->alpha) aux = aux->next;
  ins->next = aux;
  ins->prev = aux->prev;
  ins->prev->next = ins;
  ins->next->prev = ins;
}

node *create(int x, int y, node *next, node *prev, node *nextPoly,
  node *neighbor, int intersect, int entry, int visited, float alpha)
{
  node *newNode = (node*)malloc(sizeof(node));
  newNode->x = x;
  newNode->y = y;
  newNode->next = next;
  newNode->prev = prev;
  if(prev) newNode->prev->next = newNode;
  if(next) newNode->next->prev = newNode;
  newNode->nextPoly = nextPoly;
  newNode->neighbor = neighbor;
  newNode->intersect = intersect;
  newNode->entry = entry;
  newNode->visited = visited;
  newNode->alpha = alpha;
  return newNode;
}

//skip the intersection point
//if it is not a intersection point, return itself
node *next_node(node *p)
{
  node *aux=p;
  while(aux && aux->intersect) aux=aux->next;
  return aux;
}

node *last_node(node *p)
{
  node *aux=p;
  if(aux) while(aux->next) aux=aux->next;
  return aux;
}

node *first(node *p)
{
  node *aux=p;

  if (aux)
  do aux=aux->next;
  while(aux!=p && (!aux->intersect || aux->intersect && aux->visited));
  return aux;
}

//close the list to make it a loop
void circle(node *p)
{
  node *aux = last_node(p);
  aux->prev->next = p;
  p->prev = aux->prev;
  free(aux);
}

float dist(float x1, float y1, float x2, float y2)
{
  return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}


//int I(node *p1, node *p2, node *q1, node *q2,
//  float *alpha_p, float *alpha_q, int *xint, int *yint)
//this has to intersect where it is not the end of the lines
inline bool BIntersect(pt p1, pt p2, pt q1, pt q2)
{
  float  tp, tq, par;

  par = (float) ((p2.x - p1.x)*(q2.y - q1.y) -
                 (p2.y - p1.y)*(q2.x - q1.x));

  if (!par) return 0;                               /* parallel lines */
  tp = ((q1.x - p1.x)*(q2.y - q1.y) - (q1.y - p1.y)*(q2.x - q1.x))/par;
  tq = ((p2.y - p1.y)*(q1.x - p1.x) - (p2.x - p1.x)*(q1.y - p1.y))/par;

  //touching the boundary is not inside
  if(tp<=0 || tp>=1 || tq<=0 || tq>=1) return 0;

  return 1;
}

//touching boundary is also intersect
inline bool BIntersectIncludeBoundary(pt p1, pt p2, pt q1, pt q2)
{
  float  tp, tq, par;

  par = (float) ((p2.x - p1.x)*(q2.y - q1.y) -
                 (p2.y - p1.y)*(q2.x - q1.x));

  if (!par) return 0;                               /* parallel lines */
  tp = ((q1.x - p1.x)*(q2.y - q1.y) - (q1.y - p1.y)*(q2.x - q1.x))/par;
  tq = ((p2.y - p1.y)*(q1.x - p1.x) - (p2.x - p1.x)*(q1.y - p1.y))/par;

  //touching the boundary is not inside
  if(tp<0 || tp>1 || tq<0 || tq>1) return 0;

  return 1;
}

//line(p1, p2) is parallel with line(q1, q2)
inline bool parallel(pt p1, pt p2, pt q1, pt q2)
{
  float par = (float) ((p2.x - p1.x)*(q2.y - q1.y) -
                 (p2.y - p1.y)*(q2.x - q1.x));
  if(abs(par)<EPS)
      return true;
  else
      return false;
}

inline void Intersect(pt p1, pt p2, pt q1, pt q2,
        pt &pi, pt &qi)
{
    float tp, tq, par;

    par = (float) ((p2.x - p1.x)*(q2.y - q1.y) -
                   (p2.y - p1.y)*(q2.x - q1.x));

    if (!par)
        return;                               /* parallel lines */

    tp = ((q1.x - p1.x)*(q2.y - q1.y) - (q1.y - p1.y)*(q2.x - q1.x))/par;
    tq = ((p2.y - p1.y)*(q1.x - p1.x) - (p2.x - p1.x)*(q1.y - p1.y))/par;

    if(tp<0 || tp>1 || tq<0 || tq>1)
        return;

//    pi.in = true;
//    qi.in = true;
    pi.x = p1.x + tp*(p2.x - p1.x);
    pi.y = p1.y + tp*(p2.y - p1.y);
    qi.x = pi.x;
    qi.y = pi.y;

    //this can be replaced with tp and tq with care
    pi.loc = tp;// dist(p1.x, p1.y, x, y) / dist(p1.x, p1.y, p2.x, p2.y);
    qi.loc = tq;// dist(q1.x, q1.y, x, y) / dist(q1.x, q1.y, q2.x, q2.y);
}

int I(node *p1, node *p2, node *q1, node *q2,
  float *alpha_p, float *alpha_q, int *xint, int *yint)
{
  float x, y, tp, tq, t, par;

  par = (float) ((p2->x - p1->x)*(q2->y - q1->y) -
                 (p2->y - p1->y)*(q2->x - q1->x));

  if (!par) return 0;                               /* parallel lines */

  tp = ((q1->x - p1->x)*(q2->y - q1->y) - (q1->y - p1->y)*(q2->x - q1->x))/par;
  tq = ((p2->y - p1->y)*(q1->x - p1->x) - (p2->x - p1->x)*(q1->y - p1->y))/par;

  if(tp<0 || tp>1 || tq<0 || tq>1) return 0;

  x = p1->x + tp*(p2->x - p1->x);
  y = p1->y + tp*(p2->y - p1->y);

  //this can be replaced with tp and tq with care
  *alpha_p =  dist(p1->x, p1->y, x, y) / dist(p1->x, p1->y, p2->x, p2->y);
  *alpha_q = dist(q1->x, q1->y, x, y) / dist(q1->x, q1->y, q2->x, q2->y);
  *xint = (int) x;
  *yint = (int) y;

  return 1;
}

int test(node *point, node *p)
{
  node *aux, *left, i;
  int type=0;
  //(left, point) is a horizontal line segment passing "point"
  //if this line segment passes the polygon "p" for even number of times
  //the "point" is inside the polygon "p"
  left = create(0, point->y, 0, 0, 0, 0, 0, 0, 0, 0.);
  for(aux=p; aux->next; aux=aux->next)
  if(I(left, point, aux, aux->next, &i.alpha, &i.alpha, &i.x, &i.y)) type++;
  return type%2;
}

//assumption is the coordinates are all positive
inline bool testInside(pt p, trgl t)
{
    bool inside = false;
    pt left( -999, p.y);//create(0, point->y, 0, 0, 0, 0, 0, 0, 0, 0.);
    for(int i = 0; i < 3; i++)
    {
        if(BIntersect(left, p, t.p[i], t.p[(i+1)%3]))
            inside = !inside;
    }
    return inside;
}

void quit(Widget w, void *p)
{
  deleteNode(s);
  deleteNode(c);
  exit(0);
}

void redisplay(Widget w, int x, int y, void *p)
{
  node *aux, *poly;
  ClearDrawArea();

  if (aux=s)
  {
        SetColor(DRAW==1 ? WORK:S);
        while (aux->next && aux->next != s)
        {
                DrawLine(aux->x, aux->y, aux->next->x, aux->next->y);
                aux = aux->next;
        }
        if(DRAW!=1) DrawLine(aux->x, aux->y, s->x, s->y);

  }

  if (aux=c)
  {
        SetColor(DRAW==2 ? WORK:C);
        while (aux->next && aux->next != c)
        {
                DrawLine(aux->x, aux->y, aux->next->x, aux->next->y);
                aux = aux->next;
        }
        if(DRAW!=2) DrawLine(aux->x, aux->y, c->x, c->y);
  }

  if (root)
  {
        SetColor(POLY);
        for(poly = root; poly; poly = poly->nextPoly)
        {
           for(aux = poly; aux->next; aux = aux->next)
                DrawLine(aux->x, aux->y, aux->next->x, aux->next->y);
           DrawLine(aux->x, aux->y, ((node *)poly)->x, ((node *)poly)->y);
        }

        plot(s);
        plot(c);
  }
}

void add(Widget w, int which_button, int x, int y, void *data)
{
  node *newNode;
  if (!DRAW) return;

  if (which_button == 1)
  {
        newNode = (node*)malloc(sizeof(node));
        newNode->x = x;
        newNode->y = y;
        newNode->prev = 0;        /* not need to initialize with 0 after malloc ... */
        newNode->nextPoly = 0;
        newNode->neighbor = 0;
        newNode->intersect = 0;
        newNode->entry = 0;
        newNode->visited = 0;
        newNode->alpha = 0.;
        if (DRAW == 1)
        {
                newNode->next = s;
                if (s) s->prev = newNode;
                s = newNode;
        }
        else /* DRAW == 2 */
        {
                newNode->next = c;
                if (c) c->prev = newNode;
                c = newNode;
        }
        redisplay(W[3], X, Y, NULL);
  }
  else if (which_button == 3)
  {
        DRAW = DRAW==1 ? 2:0;
        redisplay(W[3], X, Y, NULL);
  }
}


void AddPoints()
{
 /*  add(NULL,1,50,50,NULL);
    add(NULL,1,200,50,NULL);
    add(NULL,1,120,150,NULL);
    add(NULL,3,10,10,NULL);

    add(NULL,1,50,200,NULL);
    add(NULL,1,200,200,NULL);
    add(NULL,1,120,50,NULL);
    add(NULL,3,50,300,NULL);*/
 /*
    add(NULL,1,100,100,NULL);
    add(NULL,1,200,100,NULL);
    add(NULL,1,100,200,NULL);
    add(NULL,3,10,10,NULL);

    add(NULL,1,150,100,NULL);
    add(NULL,1,300,150,NULL);
    add(NULL,1,50,300,NULL);
    add(NULL,3,50,300,NULL);
*/
/*
    add(NULL,1,151,52,NULL);
    add(NULL,1,102,151,NULL);
    add(NULL,1,202,151,NULL);
    add(NULL,3,10,10,NULL);

    add(NULL,1,201,102,NULL);
    add(NULL,1,152,201,NULL);
    add(NULL,1,101,103,NULL);
    add(NULL,3,50,300,NULL);
    */





    add(NULL,1,150,120,NULL);
    add(NULL,1,230,50,NULL);
    add(NULL,1,50,50,NULL);
    add(NULL,3,10,10,NULL);




    add(NULL,1,150,150,NULL);
    add(NULL,1,120,20,NULL);
    add(NULL,1,180,70,NULL);
    add(NULL,3,50,300,NULL);
}

inline bool samePrevNext(node* nd)
{
    return nd->prev->side == nd->next->side;
}
/*
bool compareByLoc(const pt &a, const pt &b)
{
    return  a.in && b.in && a.loc < b.loc;
}

bool compareByIn(const pt &a, const pt &b)
{
    return a.in && !b.in;
}
*/
//only label, not removal
//the case of prev and next are different
int RegularLabel(node* nd)
{
    int ret = -1;
    if(nd->prev->side == IN)
    {
        if(nd->next->side == ON || nd->next->side == OUT)
            ret = 0;
    }
    else if(nd->prev->side == OUT)
    {
        if(nd->next->side == ON || nd->next->side == IN)
            ret = 1;
    }
    else//if(nd->prev->side == ON)
    {
        if(nd->next->side == OUT)
            ret = 0;
        else if(nd->next->side == IN)
            ret = 1;
    }
    return ret;
}

void labelNode(node* nd)
{
    nd->entry = RegularLabel(nd);
    node* node_nei = nd->neighbor;
    if(nd->prev->side == IN && nd->next->side == IN)
    {
        if(node_nei->prev == node_nei->next)
        {
            nd->intersect = 0;
            nd->side = IN;
        }
        else
            nd->entry = 1;
    }
    else if(nd->prev->side == OUT && nd->next->side == OUT)
    {
        if(node_nei->prev == node_nei->next)
        {
            nd->intersect = 0;
            nd->side = OUT;
        }
        else
            nd->entry = 0;
    }
    else if(nd->prev->side == ON && nd->next->side == ON)
    {
        if(node_nei->prev == node_nei->next)
        {
            nd->intersect = 0;
            nd->side = IN;
        }
        else    //label opposite of neighbor
            nd->entry = (1 - RegularLabel(nd->neighbor));
    }
}

//when the loc is integer, it is the vertex,
//otherwise it is the intersection point
inline bool isVert(pt p)
{
    return (p.loc == 0 || p.loc == 1 || p.loc == 2);
}
/*
inline void GetClipped1(trgl ts, trgl tc, std::vector<pt> &clipped_vert)
{
    //c[0], c[3], s[0], s[3]
    //find the inside point in c
    pt inPt;
    for(int i = 0; i < 3; i++)
    {
        if(tc.p[i].in)
        {
            inPt = tc.p[i];
            break;
        }
    }
    //at most 4 intersection points in this case
    std::vector<pt> vert(ts.p, ts.p + 7);
    std::sort(vert.begin(), vert.end(), compareByIn);
    std::sort(vert.begin(), vert.end(), compareByLoc);
    bool prevVert = false;  //whether previous one is a vertex
    //insert the vertex from constraint triangle after an intersection point
    //in the subject triangle
    for(int i = 0; i < 6; i++)
    {
        pt cur_Pt = vert[i];
        bool curVert = isVert(cur_Pt);//whether current point is a vertex
        if(cur_Pt.in )
        {
            //if multiple adjacent points has the same loc
            //pick the last one
            if((cur_Pt.loc != vert[(i + 1) % 6].loc))
            {
                clipped_vert.push_back(cur_Pt);
                //make sure the current one is not vertex
                //because we want to insert after an intersection point
                //two inside points has to be seperated by vertex in the clipped polygon
                if(prevVert && !curVert)
                {
                    clipped_vert.push_back(inPt);
                    prevVert = false;
                }
            }
        }
        else
            break;
        prevVert = curVert;
    }
}


inline void GetClipped0(trgl t, std::vector<pt> &clipped_vert)
{
    //sort the result in s
    std::vector<pt> vert(t.p, t.p + 9);
    //the sort can be improved
    std::sort(vert.begin(), vert.end(), compareByIn);
    std::sort(vert.begin(), vert.end(), compareByLoc);
    //at most 6 vertices in the new polygon
    for(int i = 0; i < 6; i++)
    {
        if(vert[i].in)
        {
            if(vert[i].loc != vert[(i + 1) % 6].loc)
                clipped_vert.push_back(vert[i]);
        }
        else
            break;
    }
}
*/
inline void AddIntersection(trgl ts, trgl tc, pt *clipped_array, int &clipped_cnt)
{
    for(int ic = 0; ic < 3; ic++)
    {
        for(int is = 0; is < 3; is++)
        {
            pt insect_s, insect_c;
            Intersect(tc.p[ic], tc.p[(ic+1)%3], ts.p[is], ts.p[(is+1)%3 ],
                    insect_c, insect_s);

            if(insect_c.loc >= 0)
            {
                insect_c.loc += ic;
                if(clipped_cnt > 0)
                {
                    if(insect_c.loc > clipped_array[clipped_cnt - 1].loc)
                        clipped_array[clipped_cnt++] = insect_c;
                    else if(insect_c.loc < clipped_array[clipped_cnt - 1].loc)
                    {
                        clipped_array[clipped_cnt] = clipped_array[clipped_cnt - 1];
                        clipped_array[clipped_cnt - 1] = insect_c;
                        clipped_cnt++;
                    }
                    //else :insect_c.loc == clipped_vert[isect_cnt - 1].loc
                    //don't add anything
                }
                else
                {
                    clipped_array[0] = insect_c;
                    clipped_cnt++;
                }
            }
        }
    }
}

inline void swap(pt &p1, pt &p2)
{
    pt tmp;
    tmp = p1;
    p1 = p2;
    p2 = tmp;
}


struct instructSet
{
    bool doIns[12];
    instructSet()
    {
        for(int i = 0; i < 12; i++)
            doIns[i] = false;
    }
};
static instructSet stateSet[11];

void setStateInstr()
{
    stateSet[0].doIns[1] = true;

    stateSet[1].doIns[0] = true;
    stateSet[1].doIns[4] = true;

    stateSet[2].doIns[1] = true;
    stateSet[2].doIns[5] = true;

    stateSet[3].doIns[0] = true;
    stateSet[3].doIns[4] = true;
    stateSet[3].doIns[6] = true;

    stateSet[4].doIns[1] = true;
    stateSet[4].doIns[5] = true;
    stateSet[4].doIns[7] = true;

    stateSet[5].doIns[4] = true;
    stateSet[5].doIns[6] = true;
    stateSet[5].doIns[8] = true;

    stateSet[6].doIns[5] = true;
    stateSet[6].doIns[7] = true;
    stateSet[6].doIns[9] = true;

    stateSet[7].doIns[0] = true;
    stateSet[7].doIns[2] = true;
    stateSet[7].doIns[4] = true;
    stateSet[7].doIns[6] = true;

    stateSet[8].doIns[1] = true;
    stateSet[8].doIns[3] = true;
    stateSet[8].doIns[5] = true;
    stateSet[8].doIns[7] = true;

    stateSet[9].doIns[1] = true;
    stateSet[9].doIns[5] = true;
    stateSet[9].doIns[10] = true;
    stateSet[9].doIns[11] = true;

    stateSet[10].doIns[1] = true;
    stateSet[10].doIns[3] = true;
    stateSet[10].doIns[5] = true;
}


void ClipTriangle(Widget w, void *p)
{
    AddPoints();

    trgl ts, tc;
    int i = 0;
    for(node* auxs = s; auxs; auxs = auxs->next , i++)
    {
        ts.p[i].x = auxs->x;
        ts.p[i].y = auxs->y;
    }
    i = 0;
    for(node* auxc = c; auxc; auxc = auxc->next , i++)
    {
        tc.p[i].x = auxc->x;
        tc.p[i].y = auxc->y;
    }
    //mark inside or outside for the triangle vertices
    //and count the number of inside vertices
    setStateInstr();

    i = 0;
    //mark inside or outside for the triangle vertices
    //and count the number of inside vertices
    int cnt_in_s = 0, cnt_in_c = 0;
    for(i = 0; i < 3; i++)
    {
        if(tc.p[i].loc = testInside(tc.p[i], ts))
           cnt_in_c++;

        if(ts.p[i].loc = testInside(ts.p[i], tc))
            cnt_in_s++;
    }

    //make the "in" vertices in the front of the array
    int a[3] = {0, 1, 0};
    for(i = 0; i < 3; i++)
    {
        int idx = a[i];
        if(!tc.p[idx].loc && tc.p[idx + 1].loc)
            swap(tc.p[idx], tc.p[idx + 1]);
        if(!ts.p[idx].loc && ts.p[idx + 1].loc)
            swap(ts.p[idx], ts.p[idx + 1]);
    }

    bool test;
    if(1 == cnt_in_c && 1 == cnt_in_s)
    {
      //  test1 = BIntersectIncludeBoundary(ts.p[1], ts.p[2], tc.p[1], tc.p[2]);
        //if(test1)
        test = BIntersectIncludeBoundary(ts.p[1], ts.p[2], tc.p[0], tc.p[1]);
    }

    int state = -1;
    if(0 == cnt_in_c && 0 == cnt_in_s)
        state = 0;
    else if(0 == cnt_in_c && 1 == cnt_in_s)
        state = 1;
    else if(1 == cnt_in_c && 0 == cnt_in_s)
        state = 2;
    else if(0 == cnt_in_c && 2 == cnt_in_s)
        state = 3;
    else if(2 == cnt_in_c && 0 == cnt_in_s)
        state = 4;
    else if(0 == cnt_in_c && 3 == cnt_in_s)
        state = 5;
    else if(3 == cnt_in_c && 0 == cnt_in_s)
        state = 6;
    else if(1 == cnt_in_c && 2 == cnt_in_s)
        state = 7;
    else if(2 == cnt_in_c && 1 == cnt_in_s)
        state = 8;
    else if(1 == cnt_in_c && 1 == cnt_in_s && !test)
        state = 9;
    else// if(1 == cnt_in_c && 1 == cnt_in_s && !test1) and (1 == cnt_in_c && 1 == cnt_in_s && test1 && test2)
        state = 10;
    //+cs

    pt clipped_array[6];

    int clipped_cnt = 0;
    instructSet is = stateSet[state];
    if(is.doIns[0])//+sc
        AddIntersection(tc, ts, clipped_array, clipped_cnt);
    if(is.doIns[1])//+cs
        AddIntersection(ts, tc, clipped_array, clipped_cnt);
    if(is.doIns[2])//+c0-
    {
        clipped_array[clipped_cnt] = clipped_array[clipped_cnt - 1];
        clipped_array[clipped_cnt - 1] = tc.p[0];
        clipped_cnt++;
    }
    if(is.doIns[3])//+s0-
    {
        clipped_array[clipped_cnt] = clipped_array[clipped_cnt - 1];
        clipped_array[clipped_cnt - 1] = ts.p[0];
        clipped_cnt++;
    }
    if(is.doIns[4])//+s0
        clipped_array[clipped_cnt++] = ts.p[0];
    if(is.doIns[5])//+c0
        clipped_array[clipped_cnt++] = tc.p[0];
    if(is.doIns[6])//+s1
        clipped_array[clipped_cnt++] = ts.p[1];
    if(is.doIns[7])//+c1
        clipped_array[clipped_cnt++] = tc.p[1];
    if(is.doIns[8])//+s2
        clipped_array[clipped_cnt++] = ts.p[2];
    if(is.doIns[9])//+c2
        clipped_array[clipped_cnt++] = tc.p[2];
    if(is.doIns[10])//+r0
        clipped_array[clipped_cnt++] = clipped_array[0];
    if(is.doIns[11])//+r0_s0
        clipped_array[0] = ts.p[0];

    for(int i = 0; i < clipped_cnt; i++)
    {
        node* newNode = create(clipped_array[i].x, clipped_array[i].y, root, 0, 0, 0, 0, 0, 0, 0.);
        root = newNode;
    }

    redisplay(W[3], X, Y, NULL);
    CLIP=0;
}

void clip(Widget w, void *p)
{
    AddPoints();
  node *auxs, *auxc, *is, *ic;
  int xi, yi, e;
  float alpha_s, alpha_c;

  node *crt, *newNode, *old;
  int forward;

  if(DRAW || !CLIP || !s || !c) return;

  node* ss = s;
  node* cc = c;

  auxs = last_node(s);
  create(s->x, s->y, 0, auxs, 0, 0, 0, 0, 0, 0.);
  auxc = last_node(c);
  create(c->x, c->y, 0, auxc, 0, 0, 0, 0, 0, 0.);

  for(auxs = s; auxs->next; auxs = auxs->next)
  if(!auxs->intersect)
  for(auxc = c; auxc->next; auxc = auxc->next)
  if(!auxc->intersect)
  {
      node* nexts = next_node(auxs->next);
      node* nextc = next_node(auxc->next);
      if(I(auxs, nexts, auxc, nextc,
            &alpha_s, &alpha_c, &xi, &yi))
      {
            is = create(xi, yi, 0, 0, 0, 0, 1, 0, 0, alpha_s);
            ic = create(xi, yi, 0, 0, 0, 0, 1, 0, 0, alpha_c);
            is->neighbor = ic;
            ic->neighbor = is;
            //ADD-TONG-BEGIN
            is->side = ON;
            ic->side = ON;
            //ADD-TONG-END
            insert(is, auxs, next_node(auxs->next));
            insert(ic, auxc, next_node(auxc->next));

            //ADD-TONG-BEGIN
            //mark the on's
            if(alpha_s == 0)
            {
                auxs->side = ON;
            }
            else if(alpha_s == 1)
            {
                nexts->side = ON;
            }

            if(alpha_c == 0)
            {
                auxc->side = ON;
            }
            else if(alpha_c == 1)
            {
                nextc->side = ON;
            }
            //ADD-TONG-END
      }
  }

  //ADD-TONG-BEGIN
  //mark the in and out
  for(auxs = s; auxs->next; auxs = auxs->next)
  {
      if(auxs->side != ON)  //including the intersection points
      {
          if(test(auxs, c))   //inside;
              auxs->side = IN;
          else
              auxs->side = OUT;
      }
  }

  for(auxc = c; auxc->next; auxc = auxc->next)
  {
      if(auxc->side != ON)  //including the intersection points
      {
          if(test(auxc, s))   //inside;
              auxc->side = IN;
          else
              auxc->side = OUT;
      }
  }
  //ADD-TONG-END

  //MOD-TONG-FROM
  //label the entry/exit
#if TRADITIONAL
   e = test(s, c);
  if(pS) e = 1-e;
  for(auxs = s; auxs->next; auxs = auxs->next)
  if(auxs->intersect)
  {
        auxs->entry = e;
        e = 1-e;    //one entry and one exit
  }

  e=test(c, s);
  if(pC) e = 1-e;
  for(auxc = c; auxc->next; auxc = auxc->next)
  if(auxc->intersect)
  {
        auxc->entry = e;
        e = 1-e;
  }

  circle(s);
  circle(c);
#else
    //MOD-TONG-TO
  circle(s);
  circle(c);

    //label the en/ex for the intersection points
    node *heads = s;
    //skipped the head, because the head is not a intersection point
    for(auxs = heads->next; auxs != heads; auxs = auxs->next)
    {
        if(auxs->intersect)
        {
            auxc = auxs->neighbor;
            labelNode(auxs);
            labelNode(auxc);
            //if the label of (current, neighbor) is (ex, ex) or (en, en)
            /*
            if(auxs->entry == auxc->entry)
            {
                auxs->intersect = 0;
                auxc->intersect = 0;
                if(auxs->entry == 1)
                {
                    auxs->side = IN;
                    auxs->side = IN;
                }
                else//(auxs->entry == 0)
                {
                    auxs->side = OUT;
                    auxs->side = OUT;
                }
            }
            */
        }
    }

  //MOD-TONG-END
#endif

  while ((crt = first(s)) != s)
  {
        old = 0;
        for(; !crt->visited; crt = crt->neighbor)
        for(forward = crt->entry ;; )
        {
                newNode = create(crt->x, crt->y, old, 0, 0, 0, 0, 0, 0, 0.);
                old = newNode;
                crt->visited = 1;
                crt = forward ? crt->next : crt->prev;
                if(crt->intersect)
                {
                        crt->visited = 1;
                        break;
                }
        }

        old->nextPoly = root;
        root = old;
  }

  view(s);
  view(c);

  redisplay(W[3], X, Y, NULL);
  CLIP=0;
}



void set(Widget w, char *c)
{
  pS=0; pC=0;
  if(*c==65 || *c==68) pS=1;
  if(*c==65 || *c==67) pC=1;
}

void display(int argc, char **argv)
{
  if (OpenDisplay(argc, argv) == FALSE) return;

  W[0] = MakeMenu("Func");
  W[1] = MakeButton("Clip", ClipTriangle, NULL);
  W[2] = MakeButton("Quit", quit, NULL);
  W[3] = MakeDrawArea(X, Y, redisplay, NULL);
  W[4] = MakeMenuItem(W[0], "A&B ", (ButtonCB) set, (void*)"A");
  W[5] = MakeMenuItem(W[0], "A|B ", (ButtonCB) set, (void*)"B");
  W[6] = MakeMenuItem(W[0], "A/B ", (ButtonCB) set, (void*)"C");
  W[7] = MakeMenuItem(W[0], "B/A ", (ButtonCB) set, (void*)"D");

  SetWidgetPos(W[1], PLACE_RIGHT, W[0], NO_CARE, NULL);
  SetWidgetPos(W[2], PLACE_RIGHT, W[1], NO_CARE, NULL);
  SetWidgetPos(W[3], PLACE_UNDER, W[0], NO_CARE, NULL);
  SetButtonDownCB(W[3], add);

  ShowDisplay();
  GetStandardColors();
  SetBgColor(W[3], BG);
}

main(int argc, char **argv)
{
  display(argc, argv);
  MainLoop();
}
