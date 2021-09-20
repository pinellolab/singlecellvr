/*global THREE, AFRAME, Utils*/

AFRAME.registerComponent('view-indicator', {
    schema: {
        active: {type: "boolean", default: false},
    },
    init: function() {
        this.rig = document.getElementById('rig');
        this.camera = document.getElementById('player-camera');
        const formattedInfo = this.formatCameraInfo(...this.getCameraInfo());
        this.el.setAttribute('text', `value: ${formattedInfo}; color: black; side: double; width: .5;`)
    },
    tick: function() {
        if (this.data.active) {
            const formattedInfo = this.formatCameraInfo(...this.getCameraInfo());
            this.el.setAttribute('text', `value: ${formattedInfo}; color: black; side: double; width: .5;`);
        }
    },
    update: function() {
        this.el.setAttribute('visible', this.data.active);
    },
    getCameraInfo: function() {
        const posx = this.rig.getAttribute('position').x.toFixed(3);
        const posy = this.rig.getAttribute('position').y.toFixed(3);
        const posz = this.rig.getAttribute('position').z.toFixed(3);
        const rotx = this.camera.getAttribute('rotation').x.toFixed(3);
        const roty = this.camera.getAttribute('rotation').y.toFixed(3);
        const rotz = this.camera.getAttribute('rotation').z.toFixed(3);

        return [posx, posy, posz, rotx, roty, rotz];
    },
    formatCameraInfo: function(posx, posy, posz, rotx, roty, rotz) {
        return `x: ${posx}, y: ${posy}, z: ${posz}, \n
                pitch: ${rotx}, yaw: ${roty}, roll: ${rotz}`;
    }
});