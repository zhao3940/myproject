/*Wei Zhao
 * cs240
 * Assingment 3
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>


int makearg(char *input, char ***agrs);
void Hostory( int start, char listhistory[30][125]);
int issubstitution(char input[125]);
void substitution(char list[30][125], char input[125],int start);
int isalias(char input[]);
void checkalisese( char input[]);
void alias(char input[]);
void exportpath( char **agrs,int word);
void run_command(int word, char **agrs);
void piping(int word, char **agrs,int pi);
int ispiping(int word,char **agrs);

char aliases[30][5][125];
int count=0;

int main()
{
	char *p;
	int tmp, x,pi;
	int word=0;
    char input[125];
	char **agrs;// decelear 
    char list[30][125];// store history	
	pid_t childpid;
	int t=0;
	int start=0;
	FILE *fo=fopen("mshrc","r");
	
	printf("-------------wecome to msh---------------\n");
		printf("?: ");
	while(t==0)
	{
		
	    fgets(input,124,fo);
	 //  fgets(input,124,stdin); 
	 printf("inout: %s",input);
		if(input[strlen(input)-1]=='\n')
			input[strlen(input)-1]='\0';
		checkalisese( input);
		if (strcmp("!!",input)==0)// chech for "!!" command
			{	if (start%30==0)
					strcpy(input,list[29]);// if it is last iterm
				else
					strcpy(input,list[start%30-1]);
					start--;
			}
		else if (input[0] == '!'&&input[1]!='!')
		{
			tmp=0; x=1;
			
			while (input[x]!='\0')
			{
				tmp=tmp*10+input[x]-'0';
				x++; 
			}
			strcpy(input,list[tmp%30]);
			strcpy(list[tmp%30],input);
		}  
		else if(issubstitution(input)==1)// ^ ^
			{
			 substitution(list,input, start%30-1);
			 strcpy(list[start%30],input);
			}
		else
			strcpy(list[start%30],input);
		word = makearg(input, &agrs);// split input string

		pi = ispiping(word,agrs);// check if pi >=0, | command 
	   if (strcmp("exit",input)==0)// qiut from shell
			exit(0);
	/*	else if (pi>=0)
			{  
				piping( word, agrs,pi);// run piping command
				printf("?: ");
			}*/
		else if(strcmp(input,"PATH") == 0)// display path
		{
			p=getenv("PATH");
			printf("%s\n",p );
			printf("?: ");
		}
		else if(isalias(input)==1)
			{
				alias(input);
				printf("?:");
			}
		else if(strcmp(agrs[0],"export")==0)//change PATH 
			 {
				exportpath( agrs,word);
				printf("?: ");
			}
		else
		{
			childpid = fork();// creat a new process
		if(childpid<0)
		{
			printf("error\n");
		}
		else if(childpid == 0)
		{	
	
			if (strcmp ("history",input)==0)
			 {
				Hostory( start,list);
				exit(0);
			 }	
			else if (pi>=0)
			{  
				piping( word, agrs,pi);// run piping command
				printf("?: ");
			}
			else
			run_command(word,agrs);// run commands
			
		}
		else
		{  
			wait(&childpid);	
			printf("?: ");
		}
		start++;
     }	
  }
	return 0;

}

int makearg(char *input, char ***agrs)// this function is used to split string
{
	int line = 1;// at leat one word on the line
	int a = 0;
	while(input[a]==' ')
		a++;
	int m=a;
	while (input[m] !='\0')
	{
		if (input[m]==' ' && input[m+1]!=' ' && input[m+1]!='\n')// if input[] is the space character means end of a word
		line++;// ward++
		m++;
	}
	
	 agrs[0] =( char **)malloc(sizeof(char *) * line);
   	 
	int i=0;
	for( i=0; i<line; i++)  		   
	{  
		agrs[0][i] = (char *)malloc(sizeof(char) * a);		   
	}
	     int b=0;	int l =0; int z;
	for ( z=a; input[z]!='\0'; z++)
	{
		 
		if (input[z] ==' ')
		{
			if(input[z+1]!=' ' && input[z+1]!='\n')
			{
				agrs[0][l][b]='\0';// when the end of word add a '\n'
				l++;// change to other line
				b=0;
			}
		}
		else 
		{
			agrs[0][l][b]=input[z];
			b++;
		}
	}
	return line;
}
void run_command(int word,char **agrs)// run commands
{
	switch(word)// use execlp() let child process to execuate new job
			{
				case 0: break;
				case 1: 
				{
						 execlp(agrs[0], agrs[0], (char*)0);
					    	 {	
					     		printf("unknow command\n");
					    	 	exit(0);
					 		}
				}
				break;
				case 2: 
				{
					execlp(agrs[0], agrs[0], agrs[1],(char*)0);
					{	
					     		printf("unknow command\n");
					    	 	exit(0);
					 		}
				}
				break;
				case 3: 
				{
					if(execlp(agrs[0], agrs[0], agrs[1],agrs[2],(char*)0)==-1)
					{	
					     		printf("unknow command\n");
					    	 	exit(0);
					 		}
				}
				break;
				case 4: 
				{
					if(execlp(agrs[0], agrs[0], agrs[1],agrs[2],agrs[3],(char*)0)==-1)
					{	
					     		printf("unknow command\n");
					    	 	exit(0);
					 		}
				}
				case 5: 
				{
					if(execlp(agrs[0], agrs[0], agrs[1],agrs[2],agrs[3],agrs[4],(char*)0)==-1)
					{	
					     		printf("unknow command\n");
					    	 	exit(0);
					 		}
				}
				break;
				default:
                   		 {
                    		    printf("Illegal Input! More than five wards\n");
                   		 		exit(0);
                   		 }
                   		 break;
			}
}
void Hostory( int start, char listhistory[30][125])// display history, max is 30
{   int a=0;
	if (start <= 29)
	{
		while (a<=start)
		{
			printf("%s\n",listhistory[a]);
			a++;
		}
	}
	else
	{
		a= start%30 + 1; 
		int aa= a-1;
		while(a<30 )
		{
			printf("%s\n",listhistory[a]);
			a++;
		}
		a=0;
		while(a<=aa )
		{
			printf("%s\n",listhistory[aa]);
			a++;
		}
	}
	

}

int issubstitution(char input[125])// check if input has substitution command or not
{
	int a=1;
	if (input[0]=='^' && input[1]!='^')
	{
		while (input[a]!='\0')
			{
				if (input[a]=='^')
					return 1;
				a++;
			}
	}
	return 0;
}

// substitution function
void substitution(char list[30][125], char input[125],int start)
{
	int a=1; int x=0; int z=0;
	int length=0;
	char cmp[125];
	char tmplist1[125];
	char tmplist2[125];
	char newlist[125];
	while(input[a]!='^')
	{
		tmplist1[x]=input[a];
		x++;
		a++;
	}
	tmplist1[x]='\0';
//	printf(" tmp1: %s\n",tmplist1);
	x=0;
	while(input[a+1]!= '\0')
	{
		tmplist2[x]=input[a+1];
		x++;
		a++;
	}
	tmplist2[x]='\0';
//	printf(" tmplist2: %s a:%d\n",tmplist2,a);
	x=0;
	while (list[start][x]!='\0')
    {
    	if (list[start][x] == tmplist1[0])
    	{
    		z=x;
    		while(z-x < strlen(tmplist1))
    		{
    			cmp[z-x]=list[start][z];
    			z++;
    		}
			cmp[strlen(tmplist1)]='\0';
			//printf(" cmp: %s\n",cmp);
    		if (strcmp(cmp,tmplist1)==0)
    		{
    			
    			while (length < x)
    			{
    				newlist[length]=list[start][length];
    				length++;
    			}
				//printf("newlist 1: %s\n",newlist);
    			while(length<(x+strlen(tmplist2)))
    			{   int b=0;
    				newlist[length]=tmplist2[b];
    				b++;
    				length++;
    			}
				//printf("newlist 2: %s\n",newlist);
    			while (list[start][z]!='\0')
    				{
    					newlist[length]=list[start][z-1];
    					z++;
    					length++;
    				}
    				newlist[length]='\0';
    				printf("%s\n",newlist);
    			strcpy(input,newlist);
    			break;
    		}
    	}
    	x++;
    }

}

int isalias(char input[])// in my shell alias must start with no space
{
	int a=0;
	char cmp[6];
	while (a<5)
		{
			cmp[a] = input[a];
			a++;
		}
	if (strcmp(cmp,"alias")==0)
		return 1;
	return 0;
}

void alias(char input[])// store the command which has been changed 
{
	int a =5;
	int l=0; int b=0;
	while (input[a]==' ')
		a++;// skip space 
	while(input[a]!='\0')
	{  if (input[a] ==' ')
		{
			if(input[a+1]!=' ' && input[a+1]!='\n')
			{
				aliases[count][l][b]='\0';// when the end of word add a '\0'
				l++;// change to other line
				b=0;
			}
		}
		else 
		{
			aliases[count][l][b]=input[a];
			b++;
		}
		a++;
	}
	count ++;
}

// if is alias, change input to the to original command
//then run the original command 
void checkalisese( char input[])
{
	int a =0;
	int z=2; int x=0; int c=0;
	while (a < count &&a<5)
	{
		if (strcmp(input,aliases[a][0])==0)
		{
			 memset(input, '\0', 125);
			for (z=2;z<5;z++)
			{	for (x=0;aliases[a][z][x]!='\0';x++)
				{
					input[c] = aliases[a][z][x];
					c++;
				}
				input [c]=' ';
				c++;
			}
			input[c-2]='\0';
		}
	a++;
	}
	//printf("input: %s\n",input);
}

void exportpath( char **agrs,int word)
{
  char *p = getenv("PATH");
  char *new = agrs[word-1];
  strcat(p,new);
  setenv("PATH",p,1); 
  char *newp=getenv("PATH");
    printf("%s\n",newp);
}

int ispiping(int word,char **agrs)
{
	int a=0;
	while(a<word)
	{
		if( strcmp(agrs[a],"|")==0 )
		return a;
		a++;
	}
  return -1;
}

void piping(int word, char **agrs,int pi)
{
	int z=0;
	int fd[2], res;
    pid_t pid;
    res = pipe(fd);
    if(res == -1) {
        perror("Pipe error!");
        exit(1);
    }
    pid = fork();
    if(pid > 0) {
        close(fd[1]);
        dup2(fd[0], 0);
        run_command(pi,agrs);
       exit(0);
		wait(NULL);
    } else if(pid == 0) {
        dup2(fd[1], 1);
        close(fd[0]);
        char **new;
	while((pi+1)<word)
	{ strcpy(new[z],agrs[pi+1]);
	 z++;
	}
	run_command(word-pi-1,new);
	exit(0);
    } else {
        perror("Fork error!");
        exit(2);
    }
}
