#include "ai.h"
#include <QTimer>
ai::ai()
{
    QTimer *t = new QTimer();
    //setRect(1000-10,300-50,10,100);
    connect(t,SIGNAL(timeout()),this,SLOT(ai_move()));
    t->start(5);
}

void ai::ai_move()
{
   srand((unsigned)time(NULL));
   random=rand()%5-2;
   if (pos().y()>-250 && pos().y()<250)
        setPos(x(),y()+random);
    else
   {
       if(pos().y() <=-250 && random <0)
         setPos(x(),y()-random);
       else if (pos().y()<=-250 && random >0)
           setPos(x(),y()+random);
       else if (pos().y()>=250 && random >0)
           setPos(x(),y()-random);
       else if (pos().y()>=250 && random <0)
           setPos(x(),y()+random);

    }

}
