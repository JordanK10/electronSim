// Created by H.Levine with Atom Array team, Lukin group (2016).

#ifndef CAMERA_SERVER_H
#define CAMERA_SERVER_H

#include <iostream>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>



class CameraServer
{
public:
	CameraServer();
	~CameraServer();

	bool startServer();
	bool acceptConnection();


	std::vector<bool> receiveIdentifiedAtomList(int numTraps);

private:
	int serverSocket;
	int cameraSocket;

};



#endif