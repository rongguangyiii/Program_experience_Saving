#include <iostream>
#include <fstream> 
#include <vector>

int main()
{
	std::cout << "please enter your name:" << std::endl;
	std::string name;
	std::cin >> name;

	std::cout << "hello, " << name << " ! " << std::endl;

	return 0;
}