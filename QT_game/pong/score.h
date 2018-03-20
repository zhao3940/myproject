#ifndef SCORE_H
#define SCORE_H
#include <QGraphicsTextItem>
class score : public QGraphicsTextItem{
public:
    score(QGraphicsItem * parent=0);
    void increase_score11();
    void increase_score12();
    int get_score11();
    int get_score12();
private :
    int Score_value11;
    int Score_value12;

};
#endif // SCORE_H
