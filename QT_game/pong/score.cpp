#include "score.h"
#include <QFont>

score::score(QGraphicsItem *parent): QGraphicsTextItem(parent)
{
    //initialize the score to 0
    Score_value11 = 0 ;
    Score_value12 = 0 ;
    // add the text
    setPlainText(QString("Score: ")+QString::number(Score_value11)+QString(" VS ")+QString::number(Score_value12));
    setDefaultTextColor(Qt::blue);
    setFont(QFont("times",18));
}

int score::get_score11()
{
    return Score_value11;
}
int score::get_score12()
{
    return Score_value12;
}

void score::increase_score11()
{
    Score_value11++;
    setPlainText(QString("Score: ")+QString::number(Score_value11)+QString(" VS ")+QString::number(Score_value12));
    setDefaultTextColor(Qt::blue);
    setFont(QFont("times",18));
}

void score::increase_score12()
{
    Score_value12++;
    setPlainText(QString("Score: ")+QString::number(Score_value11)+QString(" VS ")+QString::number(Score_value12));
    setDefaultTextColor(Qt::blue);
    setFont(QFont("times",18));
}
