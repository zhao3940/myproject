#include "game.h"

game::game(QWidget *parent)

{
    //create sence
    Pong = new QGraphicsScene();
    //create view and set size
    QGraphicsView *view = new QGraphicsView(Pong);
    view->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    view->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    view->setFixedSize(1000,600);
    Pong->setSceneRect(0,0,1000,600);

    //set 2 player's size and postion
    player *player1=new player();
    ai *player2=new ai();
    player1->setRect(0,300-50,10,100);
    player2->setRect(1000-10,300-50,10,100);
    player1->setFlag(QGraphicsItem::ItemIsFocusable);
    player1->setFocus();
    Pong->addItem(player1);
    Pong->addItem(player2);

    //create a ball
    Ball = new ball();
    Ball->setPos(500,300);
    Pong->addItem(Ball);
    // create score bord
    sc= new score();
    Pong->addItem(sc);

    //creat help menu
    help_button *btn = new help_button(QString("Press H for Help"));
    btn->setPos(400,0);
    Pong->addItem(btn);

    //this connect() function worked first than it was not work
    // then I use KeypressEvent
    connect(btn,SIGNAL(clicked()),this,SLOT(help_if()));
    //show the view
    view->show();
}


game::~game()
{


}

void game::help_if()
{
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

