#include<iostream>  
#include<string>  
#include<stack>  
  
using namespace std;  
  
int getPriority(char ch)  
{  
    //��ȡ���ȼ�  
    if (ch == '(') return 1;  
    else if (ch == '+' || ch == '-') return 2;  
    else if (ch == '*' || ch == '/') return 3;  
    else return 4;  
}  
  
void calculate(stack<double> &mystack, char operation)  
{  
    double num1, num2, num3;  
    num2 = mystack.top();  
    mystack.pop();  
    num1 = mystack.top();  
    mystack.pop();  
    if (operation == '+') {  
        num3 = num1 + num2;  
    }  
    else if (operation == '-') {  
        num3 = num1 - num2;  
    }  
    else if (operation == '*') {  
        num3 = num1 * num2;  
    }  
    else if (operation == '/') {  
        num3 = num1 / num2;  
    }  
  
    mystack.push(num3);  
}  
  
double calculator(string str)  
{  
    //������׺���ʽ,Ĭ�������ǺϷ���  
    stack<double> mystack_number;  
    stack<char> mystack_operation;  
    int i = 0, j;  
    int size = str.size();  
    char tmp_operation;  
    string tmp_num;  
    while (i < size) {  
        if (str[i] >= '0' && str[i] <= '9') {  
            j = i;  
            while (j < size && str[j] >= '0' && str[j] <= '9') { j++; }  
            tmp_num = str.substr(i, j - i);  
            mystack_number.push(atoi(tmp_num.c_str()));  
            i = j;  
        }  
        else if (str[i] == '+' || str[i] == '-' || str[i] == '*' || str[i] == '/') {  
            if (mystack_operation.empty()) {  
                mystack_operation.push(str[i]);  
            }  
            else {  
                while (!mystack_operation.empty()) {  
                    tmp_operation = mystack_operation.top();  
                    if (getPriority(tmp_operation) >= getPriority(str[i])) {  
                        //����  
                        calculate(mystack_number, tmp_operation);  
                        mystack_operation.pop();  
                    }  
                    else break;  
                }  
                mystack_operation.push(str[i]);  
            }  
            i++;  
        }  
        else {  
            if (str[i] == '(') mystack_operation.push(str[i]);  
            else {  
                while (mystack_operation.top() != '(') {  
                    tmp_operation = mystack_operation.top();  
                    //����  
                    calculate(mystack_number, tmp_operation);  
                    mystack_operation.pop();  
                }  
                mystack_operation.pop();  
            }  
            i++;  
        }  
  
    }  
    //���������ջ�ǿգ���������Ԫ��  
    while (!mystack_operation.empty()) {  
        tmp_operation = mystack_operation.top();  
        //����  
        calculate(mystack_number, tmp_operation);  
        mystack_operation.pop();  
    }  
    return mystack_number.top();  
}  
  
int main()  
{  
    string str = "1+(2-3)*4+10/2+2*2+2+2/5";  
    cout << "��׺���ʽΪ:" << endl << str << endl;  
    double num_res = calculator(str);  
    cout << "������:" << endl << num_res << endl;  
    system("pause");  
    return 0;  
}  