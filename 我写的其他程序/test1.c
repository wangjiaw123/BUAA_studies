/*
  蔚来汽车 第三道笔试题，实现一个网络协议编解码器
*/


#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#define MAXnum 1000

int HexadecimalToInt(char s[]){
	int s_len = strlen(s);
	int out = 0;
	int i = 0;
	int tmp = 0;
	while(i<s_len){
		//printf("*i:%d,s_len:%d,s:%s\n ",i,s_len,s);
		if(s[i]>='0' && s[i]<= '9'){
			tmp = s[i] - '0';
		}
		if(s[i]>= 'A' && s[i]<= 'F'){
			tmp = s[i] - 'A'+10;
		}
		if(s[i]>= 'a' && s[i]<= 'f'){
			tmp = s[i] - 'a'+10;
		}		
		out = 16*out + tmp;
		
		i++;
		//printf("***i:%d,s_len:%d,s:%s,temp:%d***\n ",i,s_len,s,tmp);
	}
	return out;
}

int main(int argc,char **argv){
	char strinput[MAXnum];
	char str_print[MAXnum];
	int start_input = 0;
	
	/*
	while(scanf("%c",&str[start_input])!=EOF){
		start_input++;
	}
	*/
	gets(strinput);
	start_input = strlen(strinput);
	//str[start_input] = '\0';
	printf("The length of the strinput is:%d\n",start_input);
	printf("The strinput is:%s\n ",strinput);
	
	int start,end;
	start = 0;
	end = start_input;
	//53530e00420048656c6c6f2053656153756e6e6e53530b00420049666e6e702054666254766f53530b00420049666e6e702054666254766f
	while(start<end){
		//printf("<<<<<<<<<<<<<start:%d \n",start);
		int count = start;
		char temp[3];
		temp[2] = '\0';
		int len1,len2,massage_len;
		strncpy(temp,strinput+start+4,2);
		//printf("temp:%s \n",temp);
		len1 = HexadecimalToInt(temp);
		strncpy(temp,strinput+start+6,2);
		//printf("temp:%s \n",temp);
		len2 = HexadecimalToInt(temp);
		//printf("len1:%d,    len2:%d, len2<<8:%d \n",len1,len2,len2<<8);
		massage_len = len1 + (len2<<8);
		
		
		
		int total_len = 12+2*massage_len;
		//printf("start:%d, massage_len:%d ,total_len:%d \n",start,massage_len,total_len);
		if(strlen(strinput + start) < total_len){
			printf("The text is too short! Program finshed!. \n");
			return 0;
		}
		
		for(int i=0;i<2*massage_len;){
			strncpy(temp,strinput+start+12+i,2);
			//printf("//// temp = %s \n",temp);
			str_print[i/2] = (char)HexadecimalToInt(temp);
			//printf("str_print:  %s \n ",str_print);
			i+= 2;
		}
		str_print[massage_len] = '\0';
		printf("The massage is:%s \n",str_print);
		str_print[0] = '\0';
		start += total_len;
		
	}
	return 0;
}







