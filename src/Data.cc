#include "livpropa/Data.h"


namespace livpropa {


string getDataPath(string filename) {
	static string dataPath;
	if (dataPath.size())
		return concat_path(dataPath, filename);

	const char* environmentPath = getenv("LIVpropa_DATA_PATH");
	if (environmentPath) {
		if (is_directory(environmentPath)) {
			dataPath = environmentPath;
			KISS_LOG_INFO << "getDataPath: use environment variable, " << dataPath << endl;
			return concat_path(dataPath, filename);
		}
	}

	{
		string path = getInstallPrefix();
		if (path != "") {
			path += "data";
			if (is_directory(path)) {
				dataPath = path;
				KISS_LOG_INFO << "getDataPath: use install prefix, " << dataPath << endl;
				return concat_path(dataPath, filename);
			}
		}
	}

	{
		string path = executable_path() + "../data";
		if (is_directory(path)) {
			dataPath = path;
			KISS_LOG_INFO << "getDataPath: use executable path, " << dataPath << endl;
			return concat_path(dataPath, filename);
		}
	}

	dataPath = "data";
	KISS_LOG_INFO << "getDataPath: use default, " << dataPath << endl;
	return concat_path(dataPath, filename);
}


string getInstallPrefix() {
	string path = "";
	#ifdef LIVpropa_INSTALL_PREFIX
		path += LIVpropa_INSTALL_PREFIX;
	#endif

	return path;
}



} // namespace livpropa


