function save_resImage(algorithm_name)
    global ValFilename;
    fig = findobj(0,'tag','creaseg');
    ud = get(fig,'userdata');
    set(ud.txtInfo1,'string',''); %clear text: itr
    set(ud.txtInfo2,'string',''); %clear text: algo name
    imframe = getframe(get(ud.imageId,'parent'));
    imwrite(imframe.cdata,strcat('results/',ValFilename,' (',algorithm_name,').tif'));
end