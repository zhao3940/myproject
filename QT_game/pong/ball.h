#ifndef BALL_H
#define BALL_H
#include <QGraphicsRectItem>
#include <QGraphicsEllipseItem>
#include <QObject>
#include <stdlib.h>
#include <time.h>
#include <QBrush>
#include <QGraphicsTextItem>
#include <QFont>
class ball :public QObject, public QGraphicsRectItem{
    Q_OBJECT
public :
     ball();
     int random_x;
     int random_y;
public slots:
     void move();

};

#endif // BALL_H
