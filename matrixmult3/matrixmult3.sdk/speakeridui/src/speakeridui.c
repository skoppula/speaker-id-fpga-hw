/*
 * speakeridui.c: simple test application
 *
 * This application configures UART 16550 to baud rate 9600.
 * PS7 UART (Zynq) is not initialized by this application, since
 * bootrom/bsp configures it to baud rate 115200
 *
 * ------------------------------------------------
 * | UART TYPE   BAUD RATE                        |
 * ------------------------------------------------
 *   uartns550   9600
 *   uartlite    Configurable only in HW design
 *   ps7_uart    115200 (configured by bootrom/bsp)
 */

#include <stdio.h>
#include "platform.h"
#include "string.h"
#include "xil_printf.h"
#include "xgpio.h"
#include "xscugic.h"
#include "xil_exception.h"
#include "xparameters.h"
#include "xuartps.h"

#define UART_BASEADDR XPAR_PS7_UART_1_BASEADDR
#define GPIO_DEVICE_ID  XPAR_AXI_GPIO_0_DEVICE_ID	/* GPIO device that LEDs are connected to */
#define LED 0x9										/* Initial LED value - X00X */
#define LED_DELAY 40000000							/* Software delay length */
#define LED_CHANNEL 1								/* GPIO port for LEDs */

XGpio Gpio;											/* GPIO Device driver instance */

int main()
{
    init_platform();

    print("  [ZYNC] Welcome to the demo speaker ID FPGA front-end application!\n\n");

	u8 inp = 0x00;
	u32 CntrlRegister;
	CntrlRegister = XUartPs_ReadReg(UART_BASEADDR, XUARTPS_CR_OFFSET);

	XUartPs_WriteReg(UART_BASEADDR, XUARTPS_CR_OFFSET,
				  ((CntrlRegister & ~XUARTPS_CR_EN_DIS_MASK) |
				   XUARTPS_CR_TX_EN | XUARTPS_CR_RX_EN));

	volatile int Delay;
	int Status;
	Status = XGpio_Initialize(&Gpio, GPIO_DEVICE_ID);

	if (Status != XST_SUCCESS) {
			xil_printf("  [ZYNC] GPIO output to the LEDs failed!\r\n");

	} else {
		XGpio_SetDataDirection(&Gpio, LED_CHANNEL, 0x0);

		xil_printf("  [ZYNC] Press 'm' to load pre-trained DNN SID model into BRAM, 's' to start speaker verification.  Any other key to repeat menu.\n");
		do {
			XGpio_DiscreteWrite(&Gpio, LED_CHANNEL, 0);
			while (!XUartPs_IsReceiveData(UART_BASEADDR));
			inp = XUartPs_ReadReg(UART_BASEADDR, XUARTPS_FIFO_OFFSET);

			switch(inp){
			case 13:
				break;
			case 'm':
				xil_printf("  [ZYNC] Starting model load...\r\n");
				for (Delay = 0; Delay < LED_DELAY; Delay++);

				xil_printf("  [ZYNC] Waiting for model parameters over UART...\r\n");
				XGpio_DiscreteWrite(&Gpio, LED_CHANNEL, 1);
				for(int i = 0; i < 2; i++) {
					while (!XUartPs_IsReceiveData(UART_BASEADDR));
					inp = XUartPs_ReadReg(UART_BASEADDR, XUARTPS_FIFO_OFFSET);
				}
				for (Delay = 0; Delay < LED_DELAY*10; Delay++);
				xil_printf("  [ZYNC] Done receiving params and loading them to shared memory.\r\n");

				xil_printf("  [ZYNC] Press 'm' to load pre-trained DNN SID model into BRAM, 's' to start speaker verification.  Any other key to repeat menu.\n");
				break;

			case 's':
				xil_printf("  [ZYNC] Starting speaker verification...\r\n");
				for (Delay = 0; Delay < LED_DELAY; Delay++);

				xil_printf("  [ZYNC] Waiting for input samples over UART...\r\n");
				XGpio_DiscreteWrite(&Gpio, LED_CHANNEL, 8);
				for(int i = 0; i < 2; i++) {
					while (!XUartPs_IsReceiveData(UART_BASEADDR));
					inp = XUartPs_ReadReg(UART_BASEADDR, XUARTPS_FIFO_OFFSET);
				}
				xil_printf("  [ZYNC] Done receiving audio.\r\n");

				xil_printf("  [ZYNC] Computing cepstral features...\r\n");
				XGpio_DiscreteWrite(&Gpio, LED_CHANNEL, 4);
				for (Delay = 0; Delay < LED_DELAY; Delay++);
				xil_printf("  [ZYNC] Done computing cepstral features...\r\n");

				xil_printf("  [ZYNC] Evaluating frames on PL...\r\n");
				XGpio_DiscreteWrite(&Gpio, LED_CHANNEL, 2);
				for (Delay = 0; Delay < LED_DELAY*2; Delay++);
				xil_printf("  [ZYNC] Done with model evaluation.\r\n");

				xil_printf("  [ZYNC] Sending result back over UART.\n\n");
				XGpio_DiscreteWrite(&Gpio, LED_CHANNEL, 1);
				for (Delay = 0; Delay < LED_DELAY; Delay++);

				xil_printf("  [ZYNC] Press 'm' to load pre-trained DNN SID model into BRAM, 's' to start speaker verification.  Any other key to repeat menu.\n");
				break;

			default:
				xil_printf("  [ZYNC] Did not recognize character: %d. Try again! \r\n", inp);
				xil_printf("  [ZYNC] Press 'm' to load pre-trained DNN SID model into BRAM, 's' to start speaker verification.  Any other key to repeat menu.\n");
				break;
			}
		} while(1);
	}



    cleanup_platform();
    return 0;
}
