#include "livpropa/Data.h"


namespace livpropa {


std::string getDataPath(std::string filename) {
	static std::string dataPath;
	if (dataPath.size())
		return concat_path(dataPath, filename);

	const char* environmentPath = getenv("LIVPROPA_DATA_PATH");
	if (environmentPath) {
		if (is_directory(environmentPath)) {
			dataPath = environmentPath;
			KISS_LOG_INFO << "locateDataFiles: use environment variable, " << dataPath << std::endl;
			return concat_path(dataPath, filename);
		}
	}

	#ifdef LIVPROPA_INSTALL_PREFIX
	{
		std::string path = LIVPROPA_INSTALL_PREFIX "/share/livpropa";
		if (is_directory(path)) {
			dataPath = path;
			KISS_LOG_INFO << "locateDataFiles: use install prefix, " << dataPath << std::endl;
			return concat_path(dataPath, filename);
		}
	}
	#endif

	{
		std::string path = executable_path() + "../data";
		if (is_directory(path)) {
			dataPath = path;
			KISS_LOG_INFO << "locateDataFiles: use executable path, " << dataPath << std::endl;
			return concat_path(dataPath, filename);
		}
	}

	dataPath = "data";
	KISS_LOG_INFO << "locateDataFiles: use default, " << dataPath << std::endl;
	return concat_path(dataPath, filename);
}


std::string getInstallPrefix() {
	std::string path = "";
	#ifdef LIVPROPA_INSTALL_PREFIX
		path += LIVPROPA_INSTALL_PREFIX;
	#endif

	return path;
}



} // namespace livpropa


