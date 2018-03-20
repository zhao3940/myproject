#include "help_button.h"

help_button :: help_button(QString name,QGraphicsItem *parent):QGraphicsRectItem(parent)
{
    setRect(0,0,200,50);
    QBrush brush;
    brush.setStyle(Qt::SolidPattern);
    brush.setColor(Qt::darkCyan);
    setBrush(brush);
    text = new QGraphicsTextItem(name,this);
    int xPos = rect().width()/2 - text->boundingRect().width()/2;
    int yPos = rect().height()/2 - text->boundingRect().height()/2;
    text->setPos(xPos,yPos);

   // setAcceptHoverEvents(true);
}
void help_button::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
    emit clicked();
}
