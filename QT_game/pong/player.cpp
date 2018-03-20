#include "player.h"
#include <QKeyEvent>
#include <QGraphicsScene>
#include "game.h"

void player ::keyPressEvent(QKeyEvent *event)
{
    if (event->key() == Qt::Key_Up){
       if(pos().y() >-250)
        setPos(x(),y()-10);
    }
    else if (event->key() == Qt::Key_Down){
       if(pos().y()+50<300)
        setPos(x(),y()+10);
    }
    else if (event->key() == Qt::Key_H){
        QGraphicsScene *info = new QGraphicsScene();
        //create view
        QGraphicsView *view_info = new QGraphicsView(info);
        // create text
        QGraphicsTextItem *infomation1 = new QGraphicsTextItem(QString("Press Up and Down Key to control!"));
        QGraphicsTextItem *infomation2 = new QGraphicsTextItem(QString("Who got first thried scores then won!"));
        infomation2->setPos(0,10);
        info->addItem(infomation1);
        info->addItem(infomation2);
        view_info->show();
    }
}
