#ifndef AI_H
#define AI_H
#include <QGraphicsRectItem>
#include <QObject>
#include <stdlib.h>
#include <time.h>
class ai :public QObject, public QGraphicsRectItem{
    Q_OBJECT
public:
    ai();
    int random;
public slots:
    void ai_move();
};
#endif // AI_H
