#include "game.h"
#include <QApplication>

game *w;
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
     w=new game();
    //w=new game();
    //game->show();


    return a.exec();
}
