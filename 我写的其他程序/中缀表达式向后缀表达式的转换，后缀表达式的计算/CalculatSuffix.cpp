#include<iostream>  
#include<string>  
#include<stack>  
#include <iomanip>
  
using namespace std;  
  
int getPriority(char ch)  
{  
    //��ȡ���ȼ�  
    if (ch == '(') return 1;  
    else if (ch == '+' || ch == '-') return 2;  
    else if (ch == '*' || ch == '/') return 3;  
    else return 4;  
}  
  
string getPostfixExpression(string str)  
{  
    //����׺���ʽת��Ϊ��׺���ʽ  
    //Ĭ�������ǺϷ���  
    stack<char> mystack;  
    int size = str.size();  
    int i = 0,j;  
    char tmp;  
	string tmp_num;
    string res = "";  
    while (i < size) {  
        if (str[i] >= '0' && str[i] <= '9'){  

			j=i;
			//cout << "i= "<<i <<endl;
			while (j < size && str[j+1] >= '0' && str[j+1] <= '9') { j++;}  
			//cout << "j= "<<j <<endl;
			tmp_num = str.substr(i,j-i+1); 
			//cout << "tmp_num= "<<tmp_num<<endl;
			i = j;
			//cout << "i= "<<i <<endl;
			res += "#" + tmp_num + "#";
			//cout << "res " << res << endl;
			//res.push(tmp_num);  
			
            //res.push_back(str[i]);  
        }  
        else if (str[i] == '+' || str[i] == '-' || str[i] == '*' || str[i] == '/') {  
            if (mystack.empty()) {  
                mystack.push(str[i]);  
            }  
            else {  
                while (!mystack.empty()) {  
                    tmp = mystack.top();  
                    if (getPriority(tmp) >= getPriority(str[i])) {  
                        //����ջ��Ԫ��  
                        res.push_back(tmp);  
                        mystack.pop();  
                    }  
                    else break;  
                }  
                mystack.push(str[i]);  
            }  
        }  
        else {  
		if(str[i]=='('){ mystack.push(str[i]);}  
            else {  
                while (mystack.top() != '(') {  
                    tmp = mystack.top();  
                    res.push_back(tmp);  
                    mystack.pop();  
                }  
                mystack.pop();  
            }  
        }  
        i++;  
    }  
  
    //���������ջ�ǿգ���������Ԫ��  
    while (!mystack.empty()) {  
        tmp = mystack.top();  
        res.push_back(tmp);  
        mystack.pop();  
    }  
	return res;
}
  
float calculator(string str)  
{  
    //�����׺���ʽ��ֵ��Ĭ����׺���ʽ�������ֶ���һλ�ģ���0-9֮��  
    stack<float> mystack;  
    int size = str.size();  
    float num1, num2, num3;  
	int j;
	string tmp_num; 
    for (int i = 0; i < size; i++) {  
        if (str[i] == '#' && str[i+1] >= '0' && str[i+1] <= '9') {  
			j=i+1;
			while (str[j] != '#' ) { j++; } 
			tmp_num = str.substr(i+1, j-i-1);
			cout << "tmp_num= "<<tmp_num<<endl;
            mystack.push(atoi(tmp_num.c_str())); 
			i = j;	
			cout << "str[i-1:i:i+1] " << str[i-1]<< str[i] << " mystack_top: "<<mystack.top()<<endl;			
            
        } 
/*		
		else if (str[i] == '#'){
			j=i;
			while (str[j+1] != '#' ) { j++; } 
			tmp_num = str.substr(i+1, j+2);
            mystack.push(atoi(tmp_num.c_str())); 
			i = j+3;
		}
		*/
        else {  
		    
		    cout << " str[i] " <<str[i] <<endl; 
            num2 = mystack.top();  
            mystack.pop();  
            num1 = mystack.top();  
            mystack.pop();  
            if (str[i] == '+') {  
                num3 = num1 + num2;  
            }  
            else if (str[i] == '-') {  
                num3 = num1 - num2;  
            }  
            else if (str[i] == '*') {  
                num3 = num1 * num2;  
            }  
            else if (str[i] == '/') {  
                num3 = num1 / num2;  
            }  
            mystack.push(num3); 
            cout << " mystack_top: "<<mystack.top()<<endl;			
        }  
    }  
    return mystack.top();  
}  
  
   
int main()  
{  
    string str="100*3+(2-3)*2/6*3+4/2+((3*9+9/3)/2*4+9/3)*3/17"; 
    cout <<"��׺���ʽΪ:"<< endl << str << endl;  
    string res = getPostfixExpression(str);  
    cout <<"��׺���ʽΪ:"<< endl << res << endl;  
    float num_res = calculator(res);  
	cout<<setiosflags(ios::fixed)<<setprecision(5);//��Ҫͷ�ļ�#include <iomanip>
    cout <<"������:"<< endl << num_res << endl;  
    system("pause");  
    return 0;  
}  