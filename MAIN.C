//CURVEFITTING
#include<stdio.h>
#include<conio.h>
#include<math.h>
void curvefitting()
{
	int i,n;
	float a,b;
	float x[15],y[15],xy[15],x2[15];
	float sumx=0,sumy=0,sumxy=0,sumx2=0;
	clrscr();
	printf("Curve fitting.\n");
	printf("Enter the number of x and y values given.\n");
	scanf("%d",&n);
	for(i=1;i<=n;i++)
	{
		printf("Enter values of x and corresponding y.\n");
		scanf("%f%f",&x[i],&y[i]);
	}
	for(i=1;i<=n;i++)
	{
		xy[i]=x[i]*y[i];
		x2[i]=pow(x[i],2);
	}
	for(i=1;i<=n;i++)
	{
		sumx=sumx+x[i];
		sumy=sumy+y[i];
		sumxy=sumxy+xy[i];
		sumx2=sumx2+x2[i];
	}
	printf("\nMethod of least square.\n");
	printf("Linear form(y=ax+b)\n");
	printf("Using given values to form table.\n");
	printf("x\ty\txy\tx^2\n");
	for(i=1;i<=n;i++)
		printf("%.3f\t%.3f\t%.3f\t%.3f\n",x[i],y[i],xy[i],x2[i]);
	printf("\n%.3f\t%.3f\t%.3f\t%.3f\n",sumx,sumy,sumxy,sumx2);
	printf("The normal equations are:\n");
	printf("\t%.3f = %.3fa+%db\n",sumy,sumx,n);
	printf("\t%.3f = %.3f+%.3fb\n",sumxy,sumx2,sumx);
	printf("Solving these equations we get:\n");
	a=(sumxy*n-sumx*sumy)/(sumx2*n-pow(sumx,2));
	b=(sumx*sumxy-sumx2*sumy)/(pow(sumx,2)-n*sumx2);
	printf("\n\ta=%.3f and b=%.3f\n",a,b);
	printf("Therefore, Best fit = y=%.2fx+%.2f\n",a,b);
}
//HARMONIC
void harmonic()
{
	int n;
	int i;
	float x[15];
	float y[15];
	float yc[15];
	float ys[15];
	float yc2[15];
	float ys2[15];
	float sum0=0,sum1=0,sum2=0,sum3=0,sum4=0;
	float a0,a1,a2,b1,b2;
	clrscr();
	printf("Enter the number of given x and y pairs.\n");
	scanf("%d",&n);
	for(i=0;i<n;i++)
	{
		printf("Enter value of x(in degrees) & corresponding y.\n");
		scanf("%f %f",&x[i],&y[i]);
	}
	for(i=0;i<n;i++)
	{
		x[i]=(3.14159/180)*x[i];
	}
	printf("Solution:\n");
	printf("Table for computing first two harmonic:\n");
	printf("x\ty\tycosx\tysinx\tycos2x\tysin2x\n");
	for(i=0;i<n;i++)
	{
		yc[i]=cos(x[i])*y[i];
		ys[i]=sin(x[i])*y[i];
		yc2[i]=cos(2*x[i])*y[i];
		ys2[i]=sin(2*x[i])*y[i];
		printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",x[i],y[i],yc[i],ys[i],yc2[i],ys2[i]);
	}
	for(i=0;i<n;i++)
	{
		sum0=sum0+y[i];
		sum1=sum1+yc[i];
		sum2=sum2+ys[i];
		sum3=sum3+yc2[i];
		sum4=sum4+ys2[i];
	}
	printf("\n\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n\n",sum0,sum1,sum2,sum3,sum4);
	a0=(2*sum0)/n;
	a1=(2*sum1)/n;
	b1=(2*sum2)/n;
	a2=(2*sum3)/n;
	b2=(2*sum4)/n;
	printf("Hence, the values of fourier coefficents are:\n");
	printf("a0=%.3f\na1=%.3f\na2=%.3f\nb1=%.3f\nb2=%.3f\n",a0,a1,a2,b1,b2);
	printf("\nHence, the fourier series is:\n");
	printf("y = %.3f+(%.3fcosx+%.3fsinx)+(%.3fcos2x+%.3fsin2x)+.....",(a0/2),a1,b1,a2,b2);
}

//CORREL
void correlation()
{
	int n,i,g;
	float x[15],y[15];
	float X[15],Y[15];
	float XY[15];
	float X2[15],Y2[15];
	float sumX2=0,sumY2=0,sumXY=0;
	float sumx=0,sumy=0;
	float Xm,Ym;
	float r;
	float sigx,sigy;
	clrscr();
	printf("Enter the value of n.\n");
	scanf("%d",&n);
	for(i=1;i<=n;i++)
	{
		printf("Enter the values of x and correspinding y.\n");
		scanf("%f%f",&x[i],&y[i]);
	}
	for(i=1;i<=n;i++)
	{
		sumx=sumx+x[i];
		sumy=sumy+y[i];
	}
	Xm=sumx/n;
	Ym=sumy/n;
	for(i=1;i<=n;i++)
	{
		X[i]=x[i]-Xm;
		Y[i]=y[i]-Ym;
		XY[i]=X[i]*Y[i];
		X2[i]=pow(X[i],2);
		Y2[i]=pow(Y[i],2);
	}
	for(i=1;i<=n;i++)
	{
		sumX2=sumX2+X2[i];
		sumY2=sumY2+Y2[i];
		sumXY=sumXY+XY[i];
	}
	r=((sumXY)/(sqrt(sumX2)*sqrt(sumY2)));
	printf("Solution:\n");
	printf("Here, the value of n is %d.\n",n);
	printf("Hence, mean of x = Xm = %f and mean of y = Ym = %f.\n",Xm,Ym);
	printf("Table to compute coefficent of correaltion is as follows:\n");
	printf("x\ty\tX=x-Xm\tY=y-Ym\tXY\tX^2\tY^2\n");
	for(i=1;i<=n;i++)
	{
		printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",x[i],y[i],X[i],Y[i],XY[i],X2[i],Y2[i]);
	}
	printf("\n%.2f\t%.2f\t\t\t%.2f\t%.2f\t%.2f\n",sumx,sumy,sumXY,sumX2,sumY2);
	printf("\nThe value of coefficient of correlation is %.3f.\n",r);
	printf("If you wish to obtain the regression lines, Press 1 else 0.\n");
	scanf("%d",&g);
	if(g==1)
	{
		sigx=sqrt((sumX2)/n);
		sigy=sqrt((sumY2)/n);
		printf("Regression line of y on x is hence:\n");
		printf("\ty = %.3fx+%.3f\n",(r*(sigy/sigx)),(Ym-((r*(sigy/sigx))*Xm)));
		printf("Regression line of x on y is hence:\n");
		printf("\tx = %.3fy+%.3f\n",(r*(sigx/sigy)),(Xm-((r*(sigx/sigy))*Ym)));
	}
}

///LANG
void lagrange()
{
	int n,i,j,a,d;
	float x[15],y[15];
	float m=1,p=1,t=0,l,q;
	clrscr();
	printf("Lagrange's interpolation formula.\n");
	printf("Enter the number of observations.\n");
	scanf("%d",&n);
	printf("Enter the value of x and corresponding y.\n");
	for(i=0;i<n;i++)
	{
		scanf("%f%f",&x[i],&y[i]);
	}
	printf("Enter the value of x for which you wish to calculate y.\n");
	scanf("%d",&a);
	for(i=0;i<n;i++)
	{
		m=m*(a-x[i]);
	}
	l=m;
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			if(i==j)
			{
				p=p*1;
			}
			else
			{
				p=p*(x[i]-x[j]);
			}
		}
		q=l/(a-x[i]);
		t=t+((q/p)*y[i]);
		p=1;
	}
	printf("\nThe value of y for given value of x is = %.3f",t);
}


//NEWFOR
void newfor()
{
    float ax[100], ay [100], diff[100][100],z=0,l=1,x,h,p;
    int n,i,j;
    clrscr();
    printf("\nEnter the value of n:\n");
    scanf("%d",&n);

    printf("\nEnter the values in form x,y:\n");
    for (i=0;i<n;i++)
	scanf("%f %f",&ax[i],&ay[i]);
    printf("\nEnter the value of x for which the value of y is wanted: \n");
    scanf("%f",&x);


    //now making the difference table
    //calculating the 1st order of differences
    for (i=0;i<=n-1;i++)
	diff[i][1] = ay[i+1]-ay[i];

    //now calculating the second and higher order differences
    for (j=2;j<=n-1;j++)
	for(i=0;i<=(n-1)-j;i++)
	diff[i][j] = diff[i+1][j-1] - diff[i][j-1];
    for(j=1;j<=n-1;j++)
    {
    for(i=0;i<=(n-1)-j;i++)
    printf("%.2f\t",diff[i][j]);
    printf("\n");
    }
   h=ax[1]-ax[0];
    p = (x-ax[0])/h;
    printf("\n %f",p);



    for(j=1;j<n;j++)
    {

	l=(l*(p-(j-1)))/j;
	z=z+(l*diff[0][j]);

    }
    z=z+ay[0];

   printf("\n %f",z);
   getch();
}
//newton backward
void newback()
{
    float ax[100], ay[100], diff[100][100],z=0,l=1,x,h,p;
    int n,i,j,t;
    clrscr();
    printf("\nEnter the value of n:\n");
    scanf("%d",&n);

    printf("\nEnter the values in form x,y:\n");
    for (i=0;i<n;i++)
	scanf("%f %f",&ax[i],&ay[i]);
    printf("\nEnter the value of x for which the value of y is wanted: \n");
    scanf("%f",&x);


    //now making the difference table
    //calculating the 1st order of differences
    for (i=0;i<=n-1;i++)
	diff[i][1] = ay[i+1]-ay[i];

    //now calculating the second and higher order differences
    for (j=2;j<=n-1;j++)
	for(i=0;i<=(n-1)-j;i++)
	diff[i][j] = diff[i+1][j-1] - diff[i][j-1];
    for(j=1;j<=n-1;j++)
    {
    for(i=0;i<=(n-1)-j;i++)
    printf("%f\t",diff[i][j]);
    printf("\n");
    }
   h=ax[1]-ax[0];
    p = (x-ax[n-1])/h;
    printf("\n %f",p);

    t=n-2;

    for(j=1;j<n;j++)
    {

	l=(l*(p+(j-1)))/j;
	z=z+(l*diff[t][j]);
	t--;

    }
    z=z+ay[n-1];

   printf("\n %f",z);
   getch();
}
//newtons divided difference

void dd()
{
	float x[20], y[20], num, den, dd[10][10], sol=0, x1,m[20], a=1, b=0;
	int i, j, k, n;
	clrscr();
	m[0]=1;
	printf("Enter the value of x1:\n");
	scanf("%f",&x1);
	printf("Enter the number of values:\n");
	scanf("%d",&n);
	printf("Enter the vlaues of x & y:\n");
	for(i=0;i<n;i++)
	{
		scanf("%f%f",&x[i],&y[i]);
	}
	for(i=1;i<2;i++)
	{
		for(j=1;j<=n-i;j++)
		{
			num=y[j]-y[j-1];
			den=x[j]-x[j-1];
			dd[i][j]=num/den;
		}
	}
	for (i=2;i<=n-1;i++)
	{
		for(j=1;j<=n-i;j++)
		{
			num=dd[i-1][j+1] - dd[i-1][j];
			den=x[i+j-1]-x[j-1];
			dd[i][j]=num/den;
		}
	}
	printf("The values of are:\n");
	for(i=1;i<=n-1;i++)
	{
		for (j=1;j<=n-i;j++)
		{
			printf("%f\t",dd[i][j]);
		}
		printf("\n");
	}
	for (i=1;i<=n-1;i++)
	{
		for(j=1;j<i+1;j++)
		{
			m[j] = x1-x[j-1];
			for (k=j;k<j+1;k++)
			{
				a = a*m[k];
			}
		}
		b = a*dd[i][1];
		sol = sol+b;
		a=1;
	}
	sol = sol+y[0];
	printf("%f",sol);
	getch();
}

void main()
{
	int ch,ch1;
	clrscr();
	printf("WELCOME TO M3 HELPER\n");
	printf("Press 1 to Enter or 0 to Exit\n");
	scanf("%d",&ch1);
	if(ch1==1)
	{
		printf("Enter your choice");
		printf("\n1-Curve Fitting\n2-Harmonic Analysis\n3-Correlation & Regression\n4-Lagranges Interpolation\n5-Newtons Forward Interpolation\n6-Newtons Backward Interpolation\n7-Newtons Divided Difference\n");
		scanf("%d",&ch);
		switch(ch)
		{
			case 1:curvefitting();
				break;
			case 2:harmonic();
				break;
			case 3:correlation();
				break;
			case 4:lagrange();
				break;
			case 5:newfor();
				break;
			case 6:newback();
				break;
			case 7:dd();
				break;
			default:exit(0);
		}
		printf("\nThank you for using.");
	}
	else
	{
		exit(0);
	}
}