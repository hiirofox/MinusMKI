/*
  ==============================================================================

	This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
LModelAudioProcessorEditor::LModelAudioProcessorEditor(LModelAudioProcessor& p)
	: AudioProcessorEditor(&p), audioProcessor(p)
{
	// Make sure that before the constructor has finished, you've set the
	// editor's size to whatever you need it to be.
	setResizable(true, true); // 允许窗口调整大小

	setOpaque(false);  // 允许在边框外面绘制

	//setResizeLimits(64 * 11, 64 * 5, 10000, 10000); // 设置最小宽高为300x200，最大宽高为800x600
	setSize(64 * 9, 64 * 3);
	setResizeLimits(64 * 9, 64 * 3, 64 * 13, 64 * 3);

	//constrainer.setFixedAspectRatio(11.0 / 4.0);  // 设置为16:9比例
	//setConstrainer(&constrainer);  // 绑定窗口的宽高限制

	K_Pitch.setText("pitch", "");
	K_Pitch.ParamLink(audioProcessor.GetParams(), "pitch");
	addAndMakeVisible(K_Pitch);
	K_P1.setText("p1", "");
	K_P1.ParamLink(audioProcessor.GetParams(), "p1");
	addAndMakeVisible(K_P1);
	K_P2.setText("p2", "");
	K_P2.ParamLink(audioProcessor.GetParams(), "p2");
	addAndMakeVisible(K_P2);
	K_P3.setText("p3", "");
	K_P3.ParamLink(audioProcessor.GetParams(), "p3");
	addAndMakeVisible(K_P3);
	K_P4.setText("p4", "");
	K_P4.ParamLink(audioProcessor.GetParams(), "p4");
	addAndMakeVisible(K_P4);
	K_P5.setText("p5", "");
	K_P5.ParamLink(audioProcessor.GetParams(), "p5");
	addAndMakeVisible(K_P5);
	K_P6.setText("p6", "");
	K_P6.ParamLink(audioProcessor.GetParams(), "p6");
	addAndMakeVisible(K_P6);
	K_P7.setText("p7", "");
	K_P7.ParamLink(audioProcessor.GetParams(), "p7");
	addAndMakeVisible(K_P7);

	K_N1.setText("n1", "");
	K_N1.ParamLink(audioProcessor.GetParams(), "n1");
	addAndMakeVisible(K_N1);
	K_N2.setText("n2", "");
	K_N2.ParamLink(audioProcessor.GetParams(), "n2");
	addAndMakeVisible(K_N2);
	K_N3.setText("n3", "");
	K_N3.ParamLink(audioProcessor.GetParams(), "n3");
	addAndMakeVisible(K_N3);
	K_N4.setText("n4", "");
	K_N4.ParamLink(audioProcessor.GetParams(), "n4");
	addAndMakeVisible(K_N4);
	K_N5.setText("n5", "");
	K_N5.ParamLink(audioProcessor.GetParams(), "n5");
	addAndMakeVisible(K_N5);
	K_N6.setText("n6", "");
	K_N6.ParamLink(audioProcessor.GetParams(), "n6");
	addAndMakeVisible(K_N6);
	K_N7.setText("n7", "");
	K_N7.ParamLink(audioProcessor.GetParams(), "n7");
	addAndMakeVisible(K_N7);
	startTimerHz(30);

}

LModelAudioProcessorEditor::~LModelAudioProcessorEditor()
{
}

//==============================================================================
void LModelAudioProcessorEditor::paint(juce::Graphics& g)
{
	g.fillAll(juce::Colour(0x00, 0x00, 0x00));

	g.fillAll();
	g.setFont(juce::Font("FIXEDSYS", 17.0, 1));
	g.setColour(juce::Colour(0xff00ff00));;

	int w = getBounds().getWidth(), h = getBounds().getHeight();

	//g.drawText("L-MODEL Magnitudelay", juce::Rectangle<float>(32, 16, w, 16), 1);
}

void LModelAudioProcessorEditor::resized()
{
	juce::Rectangle<int> bound = getBounds();
	int x = bound.getX(), y = bound.getY(), w = bound.getWidth(), h = bound.getHeight();
	auto convXY = juce::Rectangle<int>::leftTopRightBottom;

	K_Pitch.setBounds(32 + 64 * 0, 32, 64, 64);
	K_P1.setBounds(32 + 64 * 1, 32, 64, 64);
	K_P2.setBounds(32 + 64 * 2, 32, 64, 64);
	K_P3.setBounds(32 + 64 * 3, 32, 64, 64);
	K_P4.setBounds(32 + 64 * 4, 32, 64, 64);
	K_P5.setBounds(32 + 64 * 5, 32, 64, 64);
	K_P6.setBounds(32 + 64 * 6, 32, 64, 64);
	K_P7.setBounds(32 + 64 * 7, 32, 64, 64);

	K_N1.setBounds(32 + 64 * 1, 32 + 64, 64, 64);
	K_N2.setBounds(32 + 64 * 2, 32 + 64, 64, 64);
	K_N3.setBounds(32 + 64 * 3, 32 + 64, 64, 64);
	K_N4.setBounds(32 + 64 * 4, 32 + 64, 64, 64);
	K_N5.setBounds(32 + 64 * 5, 32 + 64, 64, 64);
	K_N6.setBounds(32 + 64 * 6, 32 + 64, 64, 64);
	K_N7.setBounds(32 + 64 * 7, 32 + 64, 64, 64);
}

void LModelAudioProcessorEditor::timerCallback()
{
	repaint();
}
