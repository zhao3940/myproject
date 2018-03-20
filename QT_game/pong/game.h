#ifndef GAME_H
#define GAME_H

#include "score.h"
#include "ball.h"
#include "help_button.h"
#include "player.h"
#include "ai.h"
#include <QMainWindow>
#include <QGraphicsView>
#include <QWidget>
#include <QGraphicsScene>
#include <QGraphicsRectItem>
#include <QGraphicsEllipseItem>
#include <QKeyEvent>

class game : public QGraphicsView
{
public:
    game(QWidget *parent = 0);
    QGraphicsScene* Pong;
    score * sc;
    ball * Ball;
    ~game();

public slots:
    void help_if();
private:

};

#endif // GAME_H
