/********************************************************************
	CAPTURE_MANAGER.cpp
*********************************************************************/

#include "CAPTURE_MANAGER.h"
#include "PROJECT_INFO.h"
#include "read_write_png.h"
#include "FILE_UTILITIES.h"
#include "SCRIPT_READER.h"
#include <boost/date_time.hpp>

CAPTURE_MANAGER::CAPTURE_MANAGER()
	: last_frame_(-1)
	, auto_run_(false)
	, auto_copy_script_(false)
	, auto_capture_(false)
	, auto_video_(false)
	, auto_delete_image_(false)
	, auto_exit_(false)
	, output_common_path_(std::string())
	, output_image_path_(std::string())
	, output_video_path_(std::string())
	, output_log_path_(std::string())
	, output_script_path_(std::string())
	, output_common_file_basename_(std::string())
	, output_image_file_basename_(std::string())
	, output_video_file_basename_(std::string())
	, output_log_file_basename_(std::string())
	, output_script_file_basename_(std::string())
{
}

CAPTURE_MANAGER::~CAPTURE_MANAGER()
{
}

void CAPTURE_MANAGER::Initialize(std::string& script_abs_path)
{
	SCRIPT_READER script_reader(script_abs_path.c_str());	
	SCRIPT_BLOCK script_block = script_reader.FindBlock("SIMULATION_WORLD");

	last_frame_			= script_block.GetInteger("last_frame", -1);
	auto_run_			= script_block.GetBoolean("auto_run");	
	auto_copy_script_	= script_block.GetBoolean("auto_copy_script", true);
	auto_capture_		= script_block.GetBoolean("auto_capture");
	auto_video_			= script_block.GetBoolean("auto_video");
	auto_delete_image_	= script_block.GetBoolean("auto_delete_image", false);
	auto_exit_			= script_block.GetBoolean("auto_exit");

	// File pathe
	output_common_path_ = PROJECT_INFO::GetScriptAbsDir();
	output_image_path_  = PROJECT_INFO::GetScriptAbsDir();
	output_video_path_  = PROJECT_INFO::GetScriptAbsDir();
	output_script_path_ = PROJECT_INFO::GetScriptAbsDir();

	if(script_block.GetString("output_common_path",0))// Common path setup
	{
		output_common_path_.append(script_block.GetString("output_common_path",0));//Saving for LOG					
		output_image_path_.append(script_block.GetString("output_image_path",  "images"));//Simulation crash can leave images in the folder. So we distinguish images folder from common folder.
		output_video_path_.append(script_block.GetString("output_common_path",0));
		output_script_path_.append(script_block.GetString("output_common_path",0));
		output_common_path_.append("\\");
		output_image_path_.append("\\");
		output_video_path_.append("\\");
		output_script_path_.append("\\");
	}	
	else// Individual path setup
	{		
		output_image_path_.append(script_block.GetString("output_image_path",  "images"));
		output_video_path_.append(script_block.GetString("output_video_path",  "videos"));		
		output_script_path_.append(script_block.GetString("output_script_path", "scripts"));		
		output_image_path_.append("\\");
		output_video_path_.append("\\");
		output_script_path_.append("\\");
	}

	bool use_script_name_as_common_name = script_block.GetBoolean("use_script_name_as_common_name", false);
	// File names (excluding file extensions)	
	if(use_script_name_as_common_name)// Common file name	
	{
		output_common_file_basename_ = PROJECT_INFO::GetScriptBaseName();
		output_image_file_basename_  = output_common_file_basename_;
		output_video_file_basename_  = output_common_file_basename_;
		output_script_file_basename_ = output_common_file_basename_;
	}
	else if(script_block.GetString("output_common_file_name",0))
	{		
		output_common_file_basename_ = script_block.GetString("output_common_file_name", 0);
		output_image_file_basename_  = output_common_file_basename_;
		output_video_file_basename_  = output_common_file_basename_;
		output_script_file_basename_ = output_common_file_basename_;
	}	
	else// Individual file name
	{		
		output_image_file_basename_  = script_block.GetString("output_image_file_name",  "image");//NOTE : When image paths & names are same, running multiple simulators can cause unexpected deletion or overwriting of image files.
		output_video_file_basename_  = script_block.GetString("output_video_file_name",  "video");		
		output_script_file_basename_ = script_block.GetString("output_script_file_name", "script");		
	}

	if(script_block.GetString("simulation_name", 0))
	{
		output_image_file_basename_  = script_block.GetString("simulation_name", 0);
		output_video_file_basename_  = script_block.GetString("simulation_name", 0);
		output_script_file_basename_ = script_block.GetString("simulation_name", 0);
	}

	//make folders for image, video, log, script
	FILE_UTILITIES file_utilities;
	if(!output_image_path_.empty())
		if(!file_utilities.Directory_Exists(output_image_path_))
			file_utilities.Create_Directory(output_image_path_);		
	if(!output_video_path_.empty())
		if(!file_utilities.Directory_Exists(output_video_path_))
			file_utilities.Create_Directory(output_video_path_);
	if(!output_script_path_.empty())
		if(!file_utilities.Directory_Exists(output_script_path_))
			file_utilities.Create_Directory(output_script_path_);	
}

void CAPTURE_MANAGER::ResetCapture()
{
	Initialize(PROJECT_INFO::GetScriptAbsPath());
}

void CAPTURE_MANAGER::CopyScript()
{
	std::string script_target;		
	script_target.append(output_script_path_);
	script_target.append(output_script_file_basename_);
	script_target.append(".sim");

	if(CopyFile(PROJECT_INFO::GetScriptAbsPath().c_str(), script_target.c_str(), false))
		std::cout << "Copying script file to \"" << script_target << "\" has been succeed." << std::endl;		
	else
		std::cout << "Copying script file to \"" << script_target << "\" has been failed." << std::endl;
}

void CAPTURE_MANAGER::AutoCopyScript()
{
	if(auto_copy_script_)
		CopyScript();
}

void CAPTURE_MANAGER::CaptureImage(int current_frame, int width, int height)
{
	//NOTE : This causes images are not saved when opengl window is minimized.
	if(last_frame_ > 0 && current_frame <= last_frame_ && (current_frame == 1 || (current_frame%5) == 0))
	{
		char image_saving_pos[MAX_PATH] = {0,};
		if (current_frame == 1)
		{
			sprintf(image_saving_pos, "%s/%s_%04d.png", 
			output_image_path_.c_str(),
			output_image_file_basename_.c_str(), 
			1);
		}
		else
		{
			sprintf(image_saving_pos, "%s/%s_%04d.png", 
			output_image_path_.c_str(),
			output_image_file_basename_.c_str(), 
			(int)current_frame/5 + 1);
		}
		READ_WRITE_PNG::PngSaveImage(image_saving_pos, width, height);	
	}else if(last_frame_ < 0)
	{
		char image_saving_pos[MAX_PATH] = {0,};
		sprintf(image_saving_pos, "%s/%s_%04d.png", 
			output_image_path_.c_str(),
			output_image_file_basename_.c_str(), 
			current_frame);
		READ_WRITE_PNG::PngSaveImage(image_saving_pos, width, height);	
	}
}

void CAPTURE_MANAGER::MakeVideoAtLastFrame()
{
	if(auto_video_) 
		MakeVideo();

	if(auto_exit_) 
		exit(1);
}

void CAPTURE_MANAGER::MakeVideo()
{
	using namespace boost::posix_time;
	static std::locale loc(std::cout.getloc(), new time_facet("%Y%m%d_%H%M%S"));

	FILE_UTILITIES file_utilities;
	std::string video_encoding_command, removing_image_target;
	std::string video_file_basename = output_video_file_basename_;
	video_file_basename += "_";

	const boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::stringstream ss;
	ss.imbue(loc);
	ss << now;

	video_file_basename += ss.str();

	video_encoding_command = PROJECT_INFO::GetAppDir();
	video_encoding_command.append("ffmpeg -i \"");	
	video_encoding_command.append(output_image_path_);	
	video_encoding_command.append(output_image_file_basename_);
	video_encoding_command.append("_%04d.png\" -r 30 -g 1 -qscale 1 \"");	
	video_encoding_command.append(output_video_path_);		
	video_encoding_command.append(video_file_basename); 
	video_encoding_command.append(".avi\"");		
	file_utilities.Execute_Process(video_encoding_command);		

	std::cout<<"-end of video encoding"<<std::endl;// debug

	//removing images in images folder
	if(auto_delete_image_ == true)
	{
		char ch[100];
		for(int i=1; i <= last_frame_;i++)
		{
			removing_image_target = output_image_path_;
			removing_image_target.append(output_image_file_basename_);
			sprintf(ch, "_%04d.png", i);
			removing_image_target.append(ch);

			remove(removing_image_target.c_str());
		}
		std::cout<<"-end of deleting images"<<std::endl;// debug
	}
}