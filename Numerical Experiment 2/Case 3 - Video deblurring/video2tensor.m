function [A,fps]=video2tensor(path)
    % This function takes a video and converts it into a three-dimensional tensor
    video = VideoReader(path); 
    fps = video.FrameRate;
    numFrames = floor(video.Duration * video.FrameRate);
    frameHeight = video.Height;
    frameWidth = video.Width;
    A = zeros(frameHeight, frameWidth, numFrames); 
    frameIndex = 1;
    
    while hasFrame(video)
        A(:,:,frameIndex) = im2double(readFrame(video));              
        frameIndex = frameIndex + 1;
    end
end