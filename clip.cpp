#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include "libsx.h"
#include "node.h"
#include "math.h"
#include "stdlib.h"
#define TRADITIONAL 0

node *s=0, *c=0, *root=0;
int DRAW=1, CLIP=1, pS=1, pC=1;
Widget W[8];

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

    add(NULL,1,100,100,NULL);
    add(NULL,1,200,100,NULL);
    add(NULL,1,100,200,NULL);
    add(NULL,3,10,10,NULL);

    add(NULL,1,150,100,NULL);
    add(NULL,1,200,150,NULL);
    add(NULL,1,50,300,NULL);
    add(NULL,3,50,300,NULL);

/*
    add(NULL,1,151,52,NULL);
    add(NULL,1,102,151,NULL);
    add(NULL,1,202,151,NULL);
    add(NULL,3,10,10,NULL);

    add(NULL,1,201,102,NULL);
    add(NULL,1,152,201,NULL);
    add(NULL,1,101,103,NULL);
    add(NULL,3,50,300,NULL);*/
}

inline bool samePrevNext(node* nd)
{
    return nd->prev->side == nd->next->side;
}

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
  W[1] = MakeButton("Clip", clip, NULL);
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
