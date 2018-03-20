#ifndef HELP_BUTTON_H
#define HELP_BUTTON_H

#include <QGraphicsRectItem>
#include <QGraphicsSceneMouseEvent>
#include <QBrush>

class help_button : public QObject,public QGraphicsRectItem{
    Q_OBJECT
public:
    help_button(QString name,QGraphicsItem* parent=NULL);

    void mousePressEvent(QGraphicsSceneMouseEvent *event);
signals:
    void clicked();
private:
    QGraphicsTextItem* text;
};

#endif // HELP_BUTTON_H
