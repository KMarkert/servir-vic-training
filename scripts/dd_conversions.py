import numpy as np

def meters2dd(inPt,scale=30):
    """Function to convert meters to decimal degrees based on the approximation
    given by: https://en.wikipedia.org/wiki/Geographic_coordinate_system
    
    Args:
        inPt (list or array): A Y,X point provided in geographic coordinates
        in that order.

    Keywords:
        scale (int): Resolution of the raster value to covert into decimal 
        degrees, must be in meters.
        
    Returns:
        list: List of Y,X resolution values converted from meters to decimal degrees
        
    """

    lat = inPt[0] #get latitude value
    
    radLat = np.deg2rad(lat) #convert degree latitude to radians
    
    a = 6378137 #radius of Earth in meters
    
    ba = 0.99664719 #constant of b/a
    
    ss = np.arctan(ba*np.tan(radLat)) #calculate the reduced latitude
    
    #factor to convert meters to decimal degrees for X axis
    xfct = (np.pi/180)*a*np.cos(ss)
     
    #factor to convert meters to decimal degrees for Y axis
    yfct = (111132.92-559.82*np.cos(2*radLat)+1.175*np.cos(4*radLat)-
              0.0023*np.cos(6*radLat))
        
    #get decimal degree resolution
    ydd = scale / yfct
    xdd = scale / xfct
    
    #return list of converted resolution values
    return [ydd,xdd]

def dd2meters(inPt,scale=0.1):
    """Function to convert decimal degrees to meters based on the approximation
        given by: https://en.wikipedia.org/wiki/Geographic_coordinate_system
        
        Args:
        inPt (list or array): A Y,X point provided in geographic coordinates
        in that order.
        
        Keywords:
        scale (int): Resolution of the raster value to covert into meters,
        must be in decimal degrees.
        
        Returns:
        list: List of Y,X resolution values converted from meters to decimal degrees
        
        """
    
    lat = inPt[0] #get latitude value
    
    radLat = np.deg2rad(lat) #convert degree latitude to radians
    
    a = 6378137 #radius of Earth in meters
    
    ba = 0.99664719 #constant of b/a
    
    ss = np.arctan(ba*np.tan(radLat)) #calculate the reduced latitude
    
    #factor to convert meters to decimal degrees for X axis
    xfct = (np.pi/180)*a*np.cos(ss)
    
    #factor to convert meters to decimal degrees for Y axis
    yfct = (111132.92-559.82*np.cos(2*radLat)+1.175*np.cos(4*radLat)-
            0.0023*np.cos(6*radLat))
        
    #get meter resolution
    y_meters = scale * yfct
    x_meters = scale * xfct
    
    #return list of converted resolution values
    return [y_meters,x_meters]

def main():
    inPt = [-1,0]
    outMeters = meters2dd(inPt,100)
    outDD = dd2meters(inPt,1)
        
            
    print outMeters,outDD
                
    return

if __name__ == "__main__":
    main()


