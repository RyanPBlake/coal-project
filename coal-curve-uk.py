import pandas as pd
import math
from matplotlib import pyplot as plt
import csv

#----------------------------------------------------------------------------------------------
#This Python file contains verious methods that are used to fit two curves, hubberts curve and
#the bell curve, to the coal production of the Uk since 1700. 
#
#----------------------------------------------------------------------------------------------



#opens the file from CSV and reads the necessary data
def rawData():
    URL = "https://raw.githubusercontent.com/RyanPBlake/coal-project/main/coal-output-uk-tonnes.csv"
    df = pd.read_csv(URL)
    x = df.Year
    y = df.Coal_Output_BEIS_2020
    return x,y
#print(rawData()[1])

#rawData() is useful however there are gaps in the data. Thus a more advanced function which inerpolates the coal output
#between years is desirable. 
def xyInterpol():
    x_interpol = []
    y_interpol = []

    #start by copying rawData() exactly
    URL = "https://raw.githubusercontent.com/RyanPBlake/coal-project/main/coal-output-uk-tonnes.csv"
    df = pd.read_csv(URL)
    x = df.Year
    y = df.Coal_Output_BEIS_2020

    
    #start creating the inerpolated lists
    for i in range(len(x)):
        x_interpol.append(x[i])
        y_interpol.append(y[i])
        
        if x[i] < 1914:   #1913 is the last year which has missing data
            #create a line y = mx+c that passes through (x[i],y[i]) and (x[i+1],y[i+1])
            delta_x = x[i+1]-x[i]
            delta_y = y[i+1]-y[i]
            m = delta_y/delta_x

            #so y = y[i] + m*(x-x[i]) but x is always greater than x[i] by and integer, year.
            for year in range(1,delta_x):
                value = y[i] + m*year
                x_interpol.append(x[i] + year)
                y_interpol.append(value)

    return (x_interpol, y_interpol)
#print(xyInterpol()[1])


#open tha data but only stores up to the peak so that the first stage of comsumption can be closer studied
def dataPreMaxium():
    x = []
    y = []
    maximum = 0
    
    URL = "https://raw.githubusercontent.com/RyanPBlake/coal-project/main/coal-output-uk-tonnes.csv"
    df = pd.read_csv(URL)

    for i in range(len(df)):
        if df.Coal_Output_BEIS_2020[i] > maximum:
            x.append(df.Year[i])
            y.append(df.Coal_Output_BEIS_2020[i])
            maximum = df.Coal_Output_BEIS_2020[i]
    return x,y
#print(dataPreMaxium())

#----------------------------------------------------------------------------------------------------------------------
#now we are going to build the estimate functions, which are going be built from the x axis so that the estimate and the
#data have the same amout of data points.


#hubbert's curve is the derivative of the logistic curve f(x) = L/(1+e^(-k(x-x_0))). using this, the function first builds the
#logistic curve and then calculates the derivative using the property that f'(x) = f(x)(1-f(x)/L)*k
#hubbert's curve with coefficient = [L, k, x_0]
def hubbertsCurve(coefficients, x):
    dgdx = []
    gamma = 1/coefficients[1] # not mathematically required but allows k to be a large number instead of small
    for i in x:
        h_i = math.e**((-gamma)*(i-coefficients[2])) 
        g_i = coefficients[0]/(1+h_i)
        dgdx_i = g_i*(1-(g_i/coefficients[0]))*gamma 
        dgdx.append(dgdx_i)
    return  dgdx
#print(hubbertsCurve([990000,1,1930], rawData()[0]))


#second function is the bell curve from statistics, f(x) = Le^(-(1/(2*k^2))*(x-x_0)^2). The 1/2 exists to make k analogous
#to the variance of the function, and x_0 the mean
#bell curve with coefficient = [L, k, x_0]
def bellCurve(coefficients, x):
    y = []
    for i in x:
        h_i = -((i-coefficients[2])**2)/(2*coefficients[1]**2)
        y_i =  coefficients[0]*(math.e**(h_i))
        y.append(y_i)
    return y
#print(bellCurve([10000000000, 20, 1913], rawData()[0]))

#Cauchy curve is a third curve but is unused in the following algorithms. 
#Cauchy curve with coefficient = [L, k, x_0]
def cauchyCurve(coefficients, x):
    y = []
    for i in x:
        h_i = coefficients[1]*(i-coefficients[2])
        y_i = (coefficients[0]*coefficients[1])/(math.pi*(1+h_i**2))
        y.append(y_i)
    return y
#print(cauchyCurve([1000000,1,1913], rawData()[0]))

#---------------------------------------------------------------------------------------------------------------------------

#MSE (Mean squared error) function that calculates the error between the estimate and the data
def MSE(y,z):
    residue = 0
    N = len(y)
    for i in range(N):
        residue += (y[i] - z[i])**2
    residue = residue/N
    return(residue)
#print(MSE(rawData()[1], hubbertsCurve([25461557327, 25, 1920], rawData()[0])))


#Brute force algorithm that calculates every possible residue for a given 4d space in the coefficient domain
def bruteForce(x, y, limit):
    best_MSE = MSE(y, hubbertsCurve([1,1,1900], x))
    for L in range (1, limit):
        for k in range(1, limit):
            for x_0 in range(1, limit):
                z = hubbertsCurve([L, k, x_0], x)
                trial = MSE(y,z)
                if trial < best_MSE:
                    best_MSE = trial
                    print(best_MSE, L, k, x_0)
    return best_MSE
#print(bruteForce(rawData()[0], rawData()[1], 10))


#Function that outputs the nearby area in the coefficient domain of a given location. 
def nearArea(current_location, dt):
    foo = []
    for i in range(len(current_location)):
        temp1 = current_location[:]
        temp2 = current_location[:]
        temp1[i] = temp1[i]+dt
        temp2[i] = temp2[i]-dt
        if temp2[i] < 0 or temp2[i] == 0 : #realistically no coefficient should ever be negative or 0
            temp2[i] = dt
        #don't add the temp if it is already in foo
        if temp1 not in foo:
            foo.append(temp1)
        if temp2 not in foo:
            foo.append(temp2)
    return foo
#print(nearArea([781250000.0, 1, 1913], 1))

#a second way of finding the local area in the coefficient domain by using a hard coded set of vectors.
def nearArea2(current_location, dt):
    foo = []
    length = len(current_location)
    directions = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1], [1000, 0, 0]]
    for direction in  directions:
        temp = []
        for i in range(length):
            value = current_location[i] + direction[i]* dt
            if value < 0:
                value = current_location[i]
            temp.append(value)
        foo.append(temp)
    return foo
#print(nearArea2([1000,20,1913], 0.01))
        
#Algorithm that find the smallest residue in the neighborhood and then moves to that location and repeat
def turnBased(x, y, model, starting_location, step_length, max_turns):
    best_MSE = MSE(y, model(starting_location, x))
    current_location = starting_location

    for turn in range(max_turns):
        best_step_MSE = best_MSE
        best_step = current_location
        for step in nearArea2(current_location, step_length): # can change nearArea to nearArea2
            step_MSE = MSE(y, model(step, x))
            if step_MSE < best_step_MSE:
                best_step_MSE = step_MSE
                best_step = step
            
        current_location = best_step
        best_MSE = best_step_MSE
    return(current_location, best_MSE)
#print(turnBased(rawData()[0], rawData()[1], hubbertsCurve, [25461557327, 25.8, 1929.05], 0.01, 10000))


#for a given L, keep updating the best estimate until the only path is to update L and
#nudges L slightly in the direction to decrease the residue 
def calibrate(x, y, model, starting_location, step_length):
    best_MSE = MSE(y, model(starting_location, x))
    current_location = starting_location

    best_step = current_location
    best_step_MSE = best_MSE
    for step in nearArea(current_location, step_length):
        step_MSE = MSE(y, model(step, x))
        if step_MSE < best_MSE:
            best_step_MSE = step_MSE
            best_step = step
    while best_step[0] == current_location[0]:
        current_location = best_step
        best_MSE = best_step_MSE
        for step in nearArea(current_location, step_length):
            step_MSE = MSE(y, model(step, x))
            if step_MSE < best_step_MSE:
                best_step_MSE = step_MSE
                best_step = step
        #just in case it gets stuck in the loop which mean that the current location is the best thus return the current location
        if current_location == best_step:
                return current_location
    
    current_location = best_step    
    current_location[1] = round(current_location[1], 2)
    current_location[2] = round(current_location[2], 2)
    return current_location
#print(calibrate(rawData()[0], rawData()[1], hubbertCurve,[1000000000, 0.5, 1920], 0.01))

#bisection method that goes towards a local minimum by finding a bounded region and bisecting it
def bisectionMethod(x, y, model, starting_location, step_length, allowence):
    location = starting_location
    lower_bound = 0
    upper_bound = 0
    test = calibrate(x, y, model, location, step_length)
    
    if location[0] < test[0]:
        while location[0] < test[0]:
            storage = location[:]
            location[0] =  2*location[0]
            location[1] = test[1]
            location[2] = test[2]
            test = calibrate(x, y, model, location, step_length)
        lower_bound = storage
        upper_bound = test   
    elif location[0] > test[0]:
        while location[0] > test[0]:
            storage = location[:]
            location[0] =  location[0]/2
            location[1] = test[1]
            location[2] = test[2]
            test = calibrate(x, y, model, location, step_length)
        lower_bound = test
        upper_bound = storage

    while upper_bound[0] - lower_bound[0] > allowence:
        bisection = (upper_bound[0] + lower_bound[0])/2
        location[0] = bisection
        test = calibrate(x, y, model, location, step_length)
        if location[0] < test[0]:
            lower_bound = test
        elif location[0] > test[0]:
            upper_bound = test
        elif location[0] == test[0]:
            return test
    return test
#print()

def impericalHubbert(x, y):
    total = 0 
    for i in y:
        total += i #Sum over all coal mined 
    total += 77000000 #add proven reserves

    best_MSE = MSE(y, hubbertsCurve([total, 1, 1913], x))
    best_values = {'k': 1, 'x_0': 1913}
    for i in range(2500, 3000):
        k = i/100
        for x_0 in range(1913, 1950):
            trial_MSE = MSE(y, hubbertsCurve([total, k, x_0], x))
            if trial_MSE < best_MSE:
                best_MSE = trial_MSE
                best_values['k'] = k
                best_values['x_0'] = x_0
    return {'L': total, 'k' : best_values['k'], 'x_0' : best_values['x_0']}
#print(impericalHubbert(xyInterpol()[0], xyInterpol()[1])) # should output {'k': 0.036, 'x_0': 1922}

#same as function above but for the bell curve
def impericalBell(x, y):
    total = 0 
    for i in y:
        total += i #Sum over all coal mined 
    total += 77000000 #add proven reserves
    
    root2pi = math.sqrt(2*math.pi)
    c = 1
    a = total/(root2pi*c)
    best_MSE = MSE(y, bellCurve([a, c, 1910], x))
    best_values = {'L' : 300000000, 'k' : c, 'x_0' : 1910}
    for x_0 in range(1913, 1940):
        for k in range(350,450):
            k = k/10
            L = total/(root2pi*k)
            trial_MSE = MSE(y, bellCurve([L, k, x_0], x))
            if trial_MSE < best_MSE:
                best_MSE = trial_MSE
                best_values['L'] = L
                best_values['k'] = k
                best_values['x_0'] = x_0 
    return best_values
#print(impericalBell(xyInterpol()[0], xyInterpol()[1])) #should output {'L': 294489739.98580635, 'k': 39, 'x_0': 1922}


#function used to calculate the gradient of steepest ascent for hubbert's curve
def gradR(x, y, coefficients):
    gamma = 1/coefficients[1]
    N = len(x)
    
    drdL = 0
    drdk = 0
    drdx_0 = 0
    
    for i in range(N):
        shift = (x[i]-coefficients[2])*gamma
        h_i = math.e**(-shift)
        g_i = coefficients[0]/(1+h_i)
        dgdx_i = g_i*(1-(g_i/coefficients[0]))*gamma
    
        
        drdL += -2/(N*coefficients[0])*(y[i] - dgdx_i)*dgdx_i
        drdk += (2*(gamma**2)/N)*(y[i] - dgdx_i)*dgdx_i*(shift*(1-(2*g_i/coefficients[0]))+1)
        drdx_0 +=  (2*gamma/N)*(y[i] - dgdx_i)*dgdx_i*(1-2*g_i/coefficients[0])

        
    #drdL = -2/(N*coefficients[0]) * drdL
    #drdk = 2*(gamma**2)*drdk/N 
    #drdx_0 = 2*gamma*drdx_0/N 
    return [drdL, drdk, drdx_0]

#normalize a vector to length 1
def Norm(vector):
    normilized_vector = []
    sum_of_squares = 0
    for i in range(len(vector)):
        sum_of_squares += vector[i]**2
    y = math.sqrt(sum_of_squares)
    for i in range(len(vector)):
        normilized_vector.append(vector[i]/y)
    return normilized_vector

# gradient descent algorithm that uses the grad function to find the path of steepest descent.
def gradientDescent(x, y, starting_location, step_length, max_turns):
    location = starting_location
    for step in range(max_turns):
        foo = gradR(x, y, location)
        bar = Norm(foo)   
        for i in range(len(location)):
            location[i] = location[i] - step_length * bar[i]

            if location[1] < 0.45:
                 location[1] = 0.45
    return location 
#print(gradientDescent(rawData()[0], rawData()[1], [30000000000, 20, 1913], 1))

#------------------------------------------------------------------------------------------------------------------

#basic estimate run through turnBased() to find local minimum against raw non-interpolted data
def main():
    data = {'x' : rawData()[0], 'y' : rawData()[1]}

    fig, ax = plt.subplots()
    ax.plot(data['x'], data['y'])
    ax.plot(data['x'], hubbertsCurve([25461557327, 25, 1929.11], data['x']))
    ax.grid()

    xlabels = ax.get_xticks().astype(int)
    ax.set_xticks(xlabels)
    ax.set_xticklabels(xlabels, rotation=30)
    ax.set_xlim(1690,2030)

    ylabels = ax.get_yticks().tolist()
    ax.set_yticks(ylabels)
    ax.set_yticklabels(['{:,.0f}'.format(x) for x in ylabels], rotation=10)
    ax.set_ylim(-10000, 300000000)    

    plt.show()
#main()

#main function but ajusted for missing years creating an uneven weighting on the residue function    
def mainAjst(): 
    data = {'x' : xyInterpol()[0], 'y' : xyInterpol()[1]}
    
    hubbert_values = [29969978332.51152, 28.48, 1921.74]
    bell_values = [258092131, 44.5, 1922]
    
    hubbert_MSE = MSE(data['y'], hubbertsCurve(hubbert_values, data['x']))
    bell_MSE = MSE(data['y'], bellCurve(bell_values, data['x']))
    print('Hubbert\'s Curve has residue ', hubbert_MSE)
    print('Bell Curve has residue ', bell_MSE)
    if hubbert_MSE < bell_MSE:
        print('Hubbert\'s Curve is a better estimate')
    else:
        print('bell Curve is a better estimate')
        
    plt.plot(data['x'], data['y'])
    plt.plot(data['x'], hubbertsCurve(hubbert_values, data['x'])) #ajusted graph
    plt.plot(data['x'], bellCurve(bell_values, data['x']))
    #plt.plot(data['x'], cauchyCurve([40640000000, 0.02, 1922.1], data['x']))
    plt.grid() 
    plt.xticks(rotation = 25)
    plt.show()
#mainAjst()


#finding the best estimate for the early data set and plotting that against the whole data
def mainEarlyEstimate():
    data = {'x' : xyInterpol()[0], 'y' : xyInterpol()[1]}
    
    hubbert_MSE = MSE(dataPreMaxium()[1], hubbertsCurve([29250900000, 25, 1920.41], dataPreMaxium()[0]))
    bell_MSE = MSE(dataPreMaxium()[1], bellCurve([359935700, 53.40, 1947.54], dataPreMaxium()[0]))
    print('Before the maximum, Hubbert\'s Curve has a residue of ', hubbert_MSE)
    print('Before the maximum, the Bell Curve has a residue of ', bell_MSE)
    if hubbert_MSE < bell_MSE:
        print('Hubbert\'s Curve is a better estimate pre maximum')
    else:
        print('bell Curve is a better estimate pre maximum')

    hubbert_MSE = MSE(data['y'], hubbertsCurve([29250900000, 25, 1920.41], data['x']))
    bell_MSE = MSE(data['y'], bellCurve([359935700, 53.40, 1947.54], data['x']))
    print('Over the whole data, Hubbert\'s Curve has a residue of ', hubbert_MSE)
    print('Over the whole data, the Bell Curve has a residue of ', bell_MSE)
    if hubbert_MSE < bell_MSE:
        print('Hubbert\'s Curve is a better estimate overall')
    else:
        print('bell Curve is a better estimate overall')
     
    plt.plot(data['x'], data['y'])
    plt.plot(data['x'], hubbertsCurve([29250900000, 25, 1920.41], data['x']))
    plt.plot(data['x'], bellCurve([359935700, 53.40, 1947.54], data['x']))
    plt.grid()
    plt.show()
#mainEarlyEstimate()


def mainImp():
    data = {'x' : xyInterpol()[0], 'y' : xyInterpol()[1]}
    hubbert_values = impericalHubbert(data['x'], data['y'])
    hubbert_curve = hubbertsCurve([hubbert_values['L'], hubbert_values['k'], hubbert_values['x_0']], data['x'])
    hubbert_MSE = MSE(data['y'], hubbert_curve)
    bell_values = impericalBell(data['x'], data['y'])
    bell_curve = bellCurve([bell_values['L'], bell_values['k'], bell_values['x_0']], data['x'])
    bell_MSE = MSE(data['y'], bell_curve)
    print('Hubbert\'s curve has residue ', hubbert_MSE)
    print('Bell curve has residue ',bell_MSE)
    if hubbert_MSE < bell_MSE:
        print('Hubbert\'s curve is a better estimate')
    else:
        print('Bell curve is a better estimate')
    plt.plot(data['x'], data['y'])
    plt.plot(data['x'], hubbert_curve)
    plt.plot(data['x'], bell_curve)
    plt.grid()
    plt.show()
#mainImp()

#bisection function uses the bisection to close down on an 
def mainBisec():
    data = {'x' : xyInterpol()[0], 'y' : xyInterpol()[1]}

    hubbert_values = bisectionMethod(data['x'], data['y'], hubbertsCurve,[25000000000, 20, 1913], 0.01, 1000)
    bell_values = bisectionMethod(data['x'], data['y'], bellCurve,[25000000000, 0.5, 1913], 1, 1000000)
    print(hubbert_values)
    # should get around [32864379882.79893, 0.03, 1918.83] for hubbert's curve
    # should get around [260162352.84765625, 44, 1922] for bell curve

    hubbert_MSE = MSE(data['y'], hubbertsCurve(hubbert_values, data['x']))
    bell_MSE = MSE(data['y'], bellCurve(bell_values, data['x']))
    print(hubbert_MSE, bell_MSE)

    plt.plot(data['x'], data['y'])
    plt.plot(data['x'], hubbertsCurve(hubbert_values, data['x']))
    plt.plot(data['x'], bellCurve(bell_values, data['x']))
    
    plt.grid()
    plt.xticks(rotation = 25)
    plt.show()
#mainBisec()

def mainGD():
    data = {'x' : xyInterpol()[0], 'y' : xyInterpol()[1]}

    starting_location = [30000000000, 29.57, 1921]
    step_length = 1
    max_turns = 1000
    hubbert_values = gradientDescent(data['x'], data['y'], starting_location, step_length, max_turns)
    print(hubbert_values)

    plt.plot(data['x'], data['y'])
    plt.plot(data['x'], hubbertsCurve(hubbert_values, data['x']))

    
    plt.grid()
    plt.show()



main()
mainAjst()
mainEarlyEstimate()
mainImp()
mainBisec()
mainGD()








