#include"ball.h"
#include<QTimer>
#include<QGraphicsScene>
#include <QList>
#include "game.h"


extern game *w;

ball::ball()
{
    srand((unsigned)time(NULL));
    random_x=rand()%20-10;
    random_y=rand()%20-10;
    if (random_x==0)
        random_x++;
    if (random_y ==0)
        random_y++;
    setRect(0,0,10,10);
    QTimer *timer = new QTimer();
    connect(timer,SIGNAL(timeout()),this,SLOT(move()));
    timer->start(20);
}

void ball::move()
{
    //move ball
    setPos(x()+random_x,y()+random_y);
    if(pos().y()<=0)
        random_y*=-1;
    else if(pos().y()>=600)
        random_y*=-1;
    else if (pos().x()<0)
    {   if(w->sc->get_score12()<2)
        {   w->Ball = new ball();
            w->Ball->setPos(500,300);
            w->Pong->addItem(w->Ball);
            w->sc->increase_score12();
            delete this;

        }
        else
        {
            QGraphicsTextItem *game_over = new QGraphicsTextItem();
            game_over->setPlainText(QString("Computer Won!"));
            game_over->setDefaultTextColor(Qt::blue);
            game_over->setFont(QFont("times",18));
            game_over->setPos(450,300);
            w->Pong->addItem(game_over);
        }

    }
    else if (pos().x()>1000)
    {
        if(w->sc->get_score11()<2)
         {
                   w->Ball = new ball();
                   w->Ball->setPos(500,300);
                   w->Pong->addItem(w->Ball);
                   w->sc->increase_score11();
                   delete this;
         }
        else
        {

            QGraphicsTextItem *game_over = new QGraphicsTextItem();
            game_over->setPlainText(QString("Player Won!"));
            game_over->setDefaultTextColor(Qt::blue);
            game_over->setFont(QFont("times",18));
            game_over->setPos(450,300);
            w->Pong->addItem(game_over);
        }
    }

   QList<QGraphicsItem * > colliding_items=collidingItems();
   for (int i = 0, n= colliding_items.size();i<n;i++)
   {
       if (typeid(*colliding_items[i]) == typeid(ai)||typeid(*colliding_items[i]) == typeid(player))
           random_x*=-1;

   }

}
