// Created by H.Levine with Atom Array team, Lukin group (2016).


#include "CameraServer.h"
#include <iomanip>
#include <iostream>
#include <string>
#include <strings.h>
#include <sstream>
#include <netdb.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <vector>
#include <cstdlib>
#include <signal.h>
#include <cstring>

using namespace std;

static const char SERVER_ADDR[] = "192.168.10.2";
static const int PORT = 1234;

static const size_t buffer_size = 1024;

static const int EXIT_CODE = -1;
static const int CLOSE_CODE = -2;

static int waitingSocket;



string CLOSE_STRING("close");				// Closes connection.
string EXIT_STRING("exit");					// Ask server to shut down.
string START_ACQUISITION_STRING("start");	// Starts acquisiton of images.
string FIND_ATOMS_STRING("find");			// Find atoms in image.
string SAVE_SEQUENCE_STRING("save");		// Save images from sequence.

CameraServer::CameraServer()
{

}


CameraServer::~CameraServer()
{
	close(serverSocket);
	close(cameraSocket);	
}



bool CameraServer::startServer() {
	cout << "Starting server..." << endl;
	serverSocket = socket(AF_INET, SOCK_STREAM, 0);

	if (serverSocket < 0) {
		cout << "Unable to open socket!" << endl;
		return false;
	}


	struct sockaddr_in server_addr;

	bzero((char *) &server_addr, sizeof(server_addr));
	server_addr.sin_family = AF_INET;
	server_addr.sin_addr.s_addr = INADDR_ANY;
	server_addr.sin_port = htons(PORT);


	int yes = 1;
	setsockopt(serverSocket, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int));

	int ret = bind(serverSocket, (struct sockaddr *)&server_addr, sizeof(server_addr));
	if (ret < 0) {
		cout << "Unable to bind to socket!" << endl;
		return false;
	}

	listen(serverSocket, 1);

}


void stopHandlingSignal() {
	struct sigaction sigIntHandler;
	sigIntHandler.sa_handler = NULL;
	sigemptyset(&sigIntHandler.sa_mask);
	sigIntHandler.sa_flags = 0;

	sigaction(SIGINT, &sigIntHandler, NULL);
}


void shutdownServer(int s) {
	//shutdown(serverSocket);
	stopHandlingSignal();
}


bool CameraServer::acceptConnection() {
	struct sockaddr_in client_addr;
	socklen_t client_addr_len = sizeof(client_addr);


	waitingSocket = serverSocket;

	// Establish a sigint handler to cancel accept if necessary.
	struct sigaction sigIntHandler;
	sigIntHandler.sa_handler = shutdownServer;
	sigemptyset(&sigIntHandler.sa_mask);
	sigIntHandler.sa_flags = 0;

	sigaction(SIGINT, &sigIntHandler, NULL);

	cameraSocket = accept(serverSocket, (struct sockaddr *)&client_addr, &client_addr_len);

	if (cameraSocket < 0) {
		cout << "Error accepting connection!" << endl;
		return false;
	}
	return true;
}




vector<bool> CameraServer::receiveIdentifiedAtomList(int numTraps) {
	char buf[1024];
	bzero(buf, 1024);

	int numBytes = read(cameraSocket, buf, 1023);

	if (strlen(buf) == 0) {
		return vector<bool>();
	} else if (strcmp(buf, "exit") == 0) {
		return vector<bool>();
	}


	stringstream listString(buf);

	vector<bool> atomsPresent;
	for (int i = 0; i < numTraps; i++) {
		int present;

		listString >> present;

		if (present == 1) {
			atomsPresent.push_back(true);
		} else {
			atomsPresent.push_back(false);
		}
	}

	return atomsPresent;
}
